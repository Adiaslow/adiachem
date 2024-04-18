#!/usr/bin/env python3
# Name: Adam Murray (admurra)
# Group Members: Matteo Anicetti, Caitlyn Haley, Krisha Shetty

from __future__ import unicode_literals
import argparse
from fastaReader import FastAreader
from multiprocessing import Pool
import os
import sys
import time
from tqdm import tqdm

class tRNA:
    def __init__(self, header, sequence):
        """ Initialize a tRNA object.

        Args:
            header: header of tRNA
            sequence: sequence of tRNA

        Returns:
            None

        Raises:
            None
        """

        # Validate header
        if header == '':
            raise ValueError("Warning: Header is empty.")
        # Validate sequence
        if sequence == '':
            raise ValueError("Sequence must not be empty.")

        self.header = header.replace(" ", "")
        self.sequence = self.cleanSequence(sequence)
        self.kmerSet = self.getKmers()
        self.uniqueKmers = set()
        self.essentialKmers = set()

    def cleanSequence(self, sequence):
        """Removes non-base characters from a sequence.

        Args:
            seq: sequence

        Returns:
            sequence: sequence without alignment characters

        Raises:
            None
        """

        return "".join([c for c in sequence if c not in ['.', '_', '-']])

    def _getKmersHelper(self, args):
        """Helper function to generate k-mers from a given start index.

        Args:
            args: tuple of sequence and start index

        Returns:
            kmerSet: set of k-mers

        Raises:
            None
        """
        sequence, i = args
        return {sequence[i:j] for j in range(i+1, len(sequence)+1)}

    def _getKmersSequential(self):
        """Generate k-mers for the tRNA sequence sequentially.

        Args:
            None

        Returns:
            kmerSet: set of k-mers

        Raises:
            None
        """
        sequence = self.sequence
        kmerSet = set()
        for i in range(len(sequence)):
            kmerSet.update(self._getKmersHelper((sequence, i)))
        return kmerSet

    def getKmers(self):
        """Generate k-mers for the tRNA sequence using multiprocessing to
            speed up the process.

        Args:
            None

        Returns:
            kmerSet: set of k-mers

        Raises:
            None
        """

        sequence = self.sequence
        args = [(sequence, i) for i in range(len(sequence))]

        # Try multiprocessing to speed up the process
        try:
            with Pool() as pool:
                # Generate k-mers in parallel
                resultSets = pool.map(self._getKmersHelper, args)
            # Combine the sets from each process
            kmerSet = set().union(*resultSets)
        # Falls back on sequential processing
        except Exception as e:
            # print(f"Multiprocessing failed due to: {e}, falling back to sequential processing.")
            # Fallback to sequential processing
            kmerSet = self.getKmersSequential()

        if len(kmerSet) == 0:
            print("No k-mers found.")

        return kmerSet

    def findUniqueKmers(self, allTrnas):
        """
        Args:
            allTrnas: list of tRNA objects

        Returns:
            uniqueKmers: set of unique k-mers

        Raises:
            None
        """

        otherKmers = set().union(*(trna.kmerSet for trna in allTrnas if trna != self))
        self.uniqueKmers = self.kmerSet - otherKmers

        if len(self.uniqueKmers) == 0:
            print("No unique k-mers found.")

        return self.uniqueKmers

    def filterEssentialKmers(self):
        """ Filters our all non-essential k-mers.
        Args:
            None

        Returns:
            essentialKmers: set of essential k-mers

        Raises:
            None
        """

        # Start with all unique k-mers
        essentialKmersTemp = set(self.uniqueKmers)

        for kmer in self.uniqueKmers:
            for otherKmer in self.uniqueKmers:
                if kmer != otherKmer and kmer in otherKmer:
                    # Discard the longer k-mer
                    essentialKmersTemp.discard(otherKmer)

        self.essentialKmers = essentialKmersTemp

        if len(self.essentialKmers) == 0:
            print("No essential k-mers found.")

        return self.essentialKmers

    def countDots(self, alignedKmer):
        """Counts leading dots ('.') in the aligned k-mer.

        Args:
            alignedKmer: aligned k-mer

        Returns:
            int: number of leading dots

        Raises:
            None
        """

        initial_len = len(alignedKmer)
        stripped_len = len(alignedKmer.lstrip('.'))
        return initial_len - stripped_len

    def countBases(self, alignedKmer):
        """Counts the number of base characters (alphabets) at the beginning of the aligned k-mer.

        Args:
            alignedKmer: aligned k-mer

        Returns:
            int: number of base characters

        Raises:
            None
        """

        # Find the first occurrence of a non-base character
        for i, char in enumerate(alignedKmer):
            if not char.isalpha():
                return i

        return len(alignedKmer)  # Return the total length if all are base characters


    def alignKmers(self):
        """
        Args:
            None

        Returns:
            sortedAlignedKmers: list of aligned k-mers

        Raises:
            None
        """

        alignedKmers = []
        kmerList = list(self.essentialKmers)

        # Align k-mers
        for i in range(len(kmerList)):
            kmer = kmerList[i]
            kmerIndex = self.sequence.index(kmer)
            alignedKmer = "."*kmerIndex + kmerList[i]
            alignedKmers.append(alignedKmer)
        sortedAlignedKmers = self.sortByDotsAndSeqLength(alignedKmers)

        return sortedAlignedKmers

    def sortByDotsAndSeqLength(self, alignedKmerList):
        """Sorts aligned k-mers first by the number of leading dots, then by the sequence length.

        Args:
            alignedKmerList: list of aligned k-mers

        Returns:
            sortedAlignedKmers: list of aligned k-mers sorted by dots and sequence length

        Raises:
            None
        """

        return sorted(alignedKmerList, key=lambda kmer: (self.countDots(kmer), self.countBases(kmer)))

def main():
    """Main function.

    Args:
        None

    Returns:
        None

    Raises:
        None
    """

    # Parse arguments
    parser = argparse.ArgumentParser(description='Process a FASTA file to find unique and essential k-mers in tRNA sequences.')
    parser.add_argument("-infile", dest="fastaFile", type=str, help='Path to the input FASTA file.')
    parser.add_argument("-outfile", dest="outFile", type=str, help='Path to output file.')

    # Parse arguments
    args = parser.parse_args()
    try:
        inFile = args.fastaFile  # Use the provided FASTA file path
    except:
        raise Exception("No input file provided.")

    try:
        outFile = args.outFile  # Use the provided output file path
    except:
        raise Exception("No output file provided.")

    if inFile is None:
        inFile = "/content/drive/MyDrive/Aa-BME160/bos-tRNA-7.fa"
    if outFile is None:
        outFile = "/content/drive/MyDrive/Aa-BME160/bos-tRNA-7-essential.fa"

    # Create a FastAreader object
    reader = FastAreader(inFile)

    trnaList = []
    # Read fasta file and create tRNA objects
    for header, seq in reader.readFasta():
            trna = tRNA(header, seq)
            trnaList.append(trna)

    trnaList = sorted(trnaList, key=lambda trna: trna.header)

    # Find unique and essential k-mers for each tRNA
    for trna in trnaList:
        trna.findUniqueKmers(trnaList)
        trna.filterEssentialKmers()

    # Output unique and essential k-mers for each tRNA
    outputs = []
    # trnaList = trnaList[:5]
    for trna in tqdm(trnaList, desc="Aligning Kmers"):
        alignedAndSortedKmers = trna.alignKmers()
        # print(trna.header)
        # print(trna.sequence)
        output = trna.header + "\n"
        output += trna.sequence + "\n"
        for alignedKmer in alignedAndSortedKmers:
            output += alignedKmer + "\n"
            # print(alignedKmer)
        outputs.append(output)
        # print()

    # Write output to file
    with open(outFile, "w") as f:
        for i, trna in enumerate(trnaList):
            f.write(outputs[i])

if __name__ == "__main__":
    main()
