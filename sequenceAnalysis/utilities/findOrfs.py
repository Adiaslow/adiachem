import sys
import numpy as np
from tqdm import tqdm
from fastaReader import FastAreader

def reverseComplement(seq):
    """Generate the reverse complement of a DNA sequence.

    Args:
        seq (str): The input DNA sequence.

    Returns:
        str: The reverse complement of the input sequence.
    """
    complement = str.maketrans('ATGCN', 'TACGN')
    return seq.translate(complement)[::-1]

def getCodingFrames(seq):
    """Generate the six reading frames of a DNA sequence.

    Args:
        seq (str): The input DNA sequence.

    Returns:
        list of numpy.ndarray: A list containing numpy arrays representing the six reading frames.
    """
    revCompSeq = reverseComplement(seq)
    frames = [seq[i:] for i in range(3)] + [revCompSeq[i:] for i in range(3)]
    maxLength = max(len(frame) for frame in frames)
    paddedFrames = []
    for i, frame in enumerate(frames):
        padding_left = 'N' * (i % 3)
        padding_right = 'N' * (maxLength - len(frame) - len(padding_left))
        paddedFrame = padding_left + frame + padding_right
        paddedFrames.append(paddedFrame)
    return [np.array(list(frame)) for frame in paddedFrames]

def encodeBase(base):
    """Encode DNA bases into numerical representations.

    Args:
        base (str): The input DNA base.

    Returns:
        int: The numerical representation of the input base.
    """
    return {'N': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4}.get(base, 0)

def encodeCodingFrames(frames):
    """Encode DNA sequences into numerical representations.

    Args:
        frames (list of str): List of DNA sequences.

    Returns:
        list of numpy.ndarray: List of numpy arrays representing the numerical representations of DNA sequences.
    """
    return [np.vectorize(encodeBase)(frame) for frame in frames]

def encodeCodons(codons):
    """Encode codons into numerical representations.

    Args:
        codons (list of str): List of codons.

    Returns:
        numpy.ndarray: 2D array representing the numerical representations of codons.
    """
    return np.array([[encodeBase(base) for base in codon] for codon in codons])

def findStartsStops(encodedFrames, encodedStartCodons, encodedStopCodons):
    """Find start and stop codons in encoded frames.

    Args:
        encodedFrames (list of numpy.ndarray): List of encoded frames.
        encodedStartCodons (numpy.ndarray): Encoded start codons.
        encodedStopCodons (numpy.ndarray): Encoded stop codons.

    Returns:
        list of numpy.ndarray: List of arrays representing start and stop codons.
    """
    startsStops = [np.zeros(frame.shape, dtype=int) for frame in encodedFrames]

    for index, frame in enumerate(encodedFrames):
        for i in range(0, len(frame) - 2, 3):
            codonSlice = frame[i:i + 3]
            for startCodon in encodedStartCodons:
                if np.dot(codonSlice, startCodon) == np.dot(startCodon, startCodon):
                    startsStops[index][i] = 1
            for stopCodon in encodedStopCodons:
                if np.dot(codonSlice, stopCodon) == np.dot(stopCodon, stopCodon):
                    startsStops[index][i] = -1

    return startsStops

def extractGeneSequence(frame, startPos, stopPos, isReverse):
    """Extract gene sequence from a frame.

    Args:
        frame (str): DNA frame.
        startPos (int): Start position of the gene.
        stopPos (int): Stop position of the gene.
        isReverse (bool): Whether the frame is reversed.

    Returns:
        str: Extracted gene sequence.
    """
    if startPos is None:
        startPos = 0
    if stopPos is None:
        stopPos = len(frame)
    geneSequence = ''.join(frame[startPos:stopPos])
    geneSequence = geneSequence.replace('N', '')
    if isReverse:
        geneSequence = reverseComplement(geneSequence)
    return geneSequence

def createGene(frameIndex, startPos, stopPos, geneSequence):
    """Create a gene record.

    Args:
        frameIndex (int): Index of the frame.
        startPos (int): Start position of the gene.
        stopPos (int): Stop position of the gene.
        geneSequence (str): Gene sequence.

    Returns:
        list: Gene record containing frame index, start and stop positions, gene sequence, and gene length.
    """
    isReverse = frameIndex >= 3
    if isReverse:
        realFrame = -(frameIndex - 2)
        realStart = len(geneSequence) - stopPos + 1  # Adjust start position
        realEnd = len(geneSequence) - startPos    # Adjust end position
    else:
        realFrame = frameIndex + 1
        realStart = startPos + 1                       # Adjust start position
        realEnd = stopPos                            # No adjustment for end position
    geneLength = len(geneSequence)
    return [realFrame, abs(realStart), abs(realEnd), geneSequence, geneLength]

def getPutativeGenes(frames, startsStops, minGeneSize=100, largestORF=True):
    """Get putative genes from frames and start/stop codons.

    Args:
        frames (list of str): List of frames.
        startsStops (list of numpy.ndarray): List of arrays representing start and stop codons.
        minGeneSize (int, optional): Minimum gene size. Defaults to 100.
        largestORF (bool, optional): Whether to select the largest ORF. Defaults to True.

    Returns:
        list: List of putative genes.
    """
    allPutativeGenes = []
    for frameIndex, (frame, startsStopsArray) in enumerate(zip(frames, startsStops)):
        putativeGenes = processFrame(frameIndex, frame, startsStopsArray, minGeneSize, largestORF)
        allPutativeGenes.extend(putativeGenes)
    allPutativeGenes.sort(key=lambda x: (-x[4], x[1]))
    return allPutativeGenes

def processFrame(frameIndex, frame, startsStopsArray, minGeneSize, largestORF):
    """Process a frame to find putative genes.

    Args:
        frameIndex (int): Index of the frame.
        frame (str): DNA frame.
        startsStopsArray (numpy.ndarray): Array representing start and stop codons.
        minGeneSize (int): Minimum gene size.
        largestORF (bool): Whether to select the largest ORF.

    Returns:
        list: List of putative genes.
    """
    putativeGenes = []
    startPos = None
    for i, val in enumerate(startsStopsArray):
        if val == 1 and startPos is None:
            startPos = i
        elif val == -1:
            stopPos = i + 3
            if startPos is not None:
                geneSequence = extractGeneSequence(frame, startPos, stopPos, frameIndex >= 3)
                if len(geneSequence) >= minGeneSize:
                    geneRecord = createGene(frameIndex, startPos, stopPos, geneSequence)
                    putativeGenes.append(geneRecord)
                startPos = None
                nextStartIndex = np.where(startsStopsArray[i+1:] == 1)[0]
                if len(nextStartIndex) > 0:
                    nextStartIndex = nextStartIndex[0] + i + 1
                    if (nextStartIndex - i) % 3 == 0:
                        startPos = nextStartIndex
            else:
                geneSequence = extractGeneSequence(frame, 0, stopPos, frameIndex >= 3)
                if len(geneSequence) >= minGeneSize:
                    geneRecord = createGene(frameIndex, 0, stopPos, geneSequence)
                    putativeGenes.append(geneRecord)
    if startPos is not None:
        stopPos = len(frame)
        geneSequence = extractGeneSequence(frame, startPos, stopPos, frameIndex >= 3)
        if len(geneSequence) >= minGeneSize:
            geneRecord = createGene(frameIndex, startPos, stopPos, geneSequence)
            putativeGenes.append(geneRecord)
    if not largestORF:
        return putativeGenes
    else:
        return [max(putativeGenes, key=lambda x: x[4])] if putativeGenes else []

def findOrfsInSequence(sequence, startCodons, stopCodons):
    """Find ORFs in a DNA sequence.

    Args:
        sequence (str): DNA sequence.
        startCodons (list of str): Start codons.
        stopCodons (list of str): Stop codons.

    Returns:
        list of numpy.ndarray: List of arrays representing start and stop codons.
    """
    frames = getCodingFrames(sequence)
    encodedFrames = encodeCodingFrames(frames)
    encodedStartCodons = encodeCodons(startCodons)
    encodedStopCodons = encodeCodons(stopCodons)
    startsStops = findStartsStops(encodedFrames, encodedStartCodons, encodedStopCodons)

    return startsStops

def formatGeneForPrinting(geneRecord):
    """Format a gene record for printing.

    Args:
        geneRecord (list): Gene record containing frame, start, end, sequence, and length.

    Returns:
        str: Formatted gene string.
    """
    frame, start, end, sequence, length = geneRecord
    formattedGene = f"{frame:+d} {start:>5d}..{end:>5d} {length:>5d}" # {sequence}"
    return formattedGene

def findOrfsInGenome(sequence, startCodons=["ATG"], stopCodons=["TAG", "TAA", "TGA"], minGeneSize=100, largestORF=True, chunk_size=10000):
    """Find ORFs in a genome.

    Args:
        sequence (str): Genome sequence.
        startCodons (list of str, optional): Start codons. Defaults to ["ATG"].
        stopCodons (list of str, optional): Stop codons. Defaults to ["TAG", "TAA", "TGA"].
        minGeneSize (int, optional): Minimum gene size. Defaults to 100.
        largestORF (bool, optional): Whether to select the largest ORF. Defaults to True.
        chunk_size (int, optional): Size of each sequence chunk to process. Defaults to 1000.

    Returns:
        list: List of found ORFs.
    """
    # Split the sequence into chunks
    chunks = [sequence[i:i+chunk_size] for i in range(0, len(sequence), chunk_size)]

    allGenes = []
    for chunk in chunks:
        forwardFrames = getCodingFrames(chunk)
        reverseSequence = reverseComplement(chunk)
        reverseFrames = getCodingFrames(reverseSequence)

        encodedStartCodons = encodeCodons(startCodons)
        encodedStopCodons = encodeCodons(stopCodons)

        forwardStartsStops = findStartsStops(encodeCodingFrames(forwardFrames), encodedStartCodons, encodedStopCodons)
        reverseStartsStops = findStartsStops(encodeCodingFrames(reverseFrames), encodedStartCodons, encodedStopCodons)

        forwardGenes = getPutativeGenes(forwardFrames, forwardStartsStops, minGeneSize, largestORF)
        reverseGenes = getPutativeGenes(reverseFrames, reverseStartsStops, minGeneSize, largestORF)

        allGenes.extend(forwardGenes + reverseGenes)

    return allGenes


def formatGenesForPrinting(genesWithHeaders):
    """Format genes with headers for printing.

    Args:
        genesWithHeaders (list of tuple): List of tuples containing header and genes.

    Returns:
        str: Formatted genes for printing.
    """
    formattedGenes = []
    for header, foundOrfs in genesWithHeaders:
        formattedGeneGroup = []
        formattedGeneGroup.append(header + "\n")

        # Merge forward and reverse genes
        allGenes = foundOrfs[::2] + foundOrfs[1::2]

        # Sort genes by decreasing ORF size
        allGenes.sort(key=lambda x: (-x[4], x[1]))

        for gene in allGenes:
            formattedGene = formatGeneForPrinting(gene)
            formattedGeneGroup.append(formattedGene + "\n")
        formattedGenes.append(formattedGeneGroup)
    return ''.join([''.join(gene) for gene in formattedGenes])

def main():
    # reader = FastAreader('/content/drive/MyDrive/Aa-BME160/lab5test.fa')
    reader = FastAreader('../tass2.fa')


    sequencesWithHeaders = [(header, seq) for header, seq in reader.readFasta()]

    loggingEnabled = True

    allFoundOrfsWithHeaders = []
    for header, seq in tqdm(sequencesWithHeaders):
        foundOrfs = findOrfsInGenome(seq, ["ATG"], ["TAG", "TAA", "TGA"], minGeneSize=100, largestORF=True)
        allFoundOrfsWithHeaders.append((header, foundOrfs))

    resultForPrinting = formatGenesForPrinting(allFoundOrfsWithHeaders)

    with open('tass2ORFdata-ATG-100.txt', 'w') as f:
        print()
        print(resultForPrinting)

if __name__ == "__main__":
    main()
