import sqlite3
import uuid

class NucParams:
    """A class for analyzing nucleotide and codon parameters from sequences.

    Attributes:
        rnaCodonTable: A mapping of RNA codons to their corresponding amino acids.
        dnaCodonTable: A mapping of DNA codons to their corresponding amino acids,
            derived from the RNA codon table.
        genomeFeatures: Stores aggregated data on nucleotide and codon composition,
            individual sequence information, and total counts.
    """

    rnaCodonTable = {
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }

    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__(self, inString=''):
        """Initializes the NucParams instance with optional sequence data.

        Args:
            inString: An optional initial nucleotide sequence to analyze.

        Returns:
            None

        Raises:
            None
        """

        # Initializes the root of the genome features dictionary
        self.genomeFeatures = {
            "Total Amino Acid Composition": {aa:0 \
                for aa in set(self.rnaCodonTable.values())},
            "Total Codon Composition": {codon:0 \
                for codon in self.rnaCodonTable.keys()},
            "Total Codon Count": 0,
            "Total Nucleic Acid Composition": {nuc:0 \
                for nuc in ['A', 'C', 'G', 'T', 'U', 'N']},
            "Total Nucleic Acid Count": 0,
            "Individual Sequences": {}
        }
        # Adds and characterizes the first sequence
        if inString:
            self.addSequence(inString)


    def _validateNuc(self, nuc):
        """Validates if a nucleotide is part of the nucleic acid composition.

        Args:
            nuc: The nucleotide to validate.

        Returns:
            True if the nucleotide is valid, False otherwise.

        Raises:
            None
        """

        return nuc.upper() in self.genomeFeatures["Total Nucleic Acid Composition"].keys()


    def _validateCodon(self, codon):
        """Validates if a codon is part of the codon composition.

        Args:
            codon: The codon to validate.

        Returns:
            True if the codon is valid, False otherwise.

        Raises:
            None
        """

        return codon in self.genomeFeatures["Total Codon Composition"].keys()


    def _getAAComp(self, seq):
        """Calculates the amino acid composition of a given sequence.

        Args:
            seq: A nucleotide sequence.

        Returns:
            A dictionary mapping each amino acid to its count in the sequence.

        Raises:
            None
        """

        comp = {aa:0 for aa in self.rnaCodonTable.values()}
        count = 0
        for codon in [seq[i:i+3] for i in range(0, len(seq), 3)]:
            codon = codon.replace('T','U')
            if self._validateCodon(codon):
                comp[self.rnaCodonTable[codon]]+=1
            else:
                print("Invalid codon.")

        return comp


    def _getCodonComp(self, seq):
        """Calculates the codon composition of a given sequence.

        Args:
            seq: A nucleotide sequence.

        Returns:
            A dictionary mapping each codon to its count in the sequence.

        Raises:
            None
        """

        comp = {codon:0 for codon in self.rnaCodonTable.keys()}
        count = 0
        for codon in [seq[i:i+3] for i in range(0, len(seq), 3)]:
            codon = codon.replace('T','U')
            # print(codon)
            if self._validateCodon(codon):
                comp[codon]+=1
            else:
                print("Invalid codon.")

        return comp


    def _getNucComp(self, seq):
        """Calculates the nucleotide composition of a given sequence.

        Args:
            seq: A nucleotide sequence.

        Returns:
            A dictionary mapping each nucleotide to its count in the sequence.

        Raises:
            None
        """

        comp = {nuc:0 for nuc in ['A', 'C', 'G', 'T', 'U']}
        count = 0
        for nuc in seq:
            if self._validateNuc(nuc):
                comp[nuc]+=1
            else:
                print("Invalid nuc.")
        return comp


    def generateSequenceId(self):
        """Generates a unique identifier for a sequence.

        Args:
            None

        Returns:
            A unique sequence identifier.

        Raises:
            None
        """

        sequenceId = str(uuid.uuid4())
        return sequenceId


    def addSequence(self, inSeq):
        """Processes and adds a nucleotide sequence to the genome features.

        This method updates the total amino acid composition, total
            codon composition, total nucleic acid composition, and individual
            sequences data with the new sequence.

        Args:
            inSeq: A nucleotide sequence to be analyzed and added.

        Returns:
            None

        Raises:
            None
        """

        aaComp = self._getAAComp(inSeq)
        codonComp = self._getCodonComp(inSeq)
        nucComp = self._getNucComp(inSeq)

        # Add individual sequence data to "Individual Sequences"
        self.genomeFeatures["Individual Sequences"][self.generateSequenceId()] = {
            "Input Sequence": inSeq,
            "Amino Acid Composition": aaComp,
            "Codon Composition": codonComp,
            "Codon Count": sum(codonComp.values()),
            "Nucleic Acid Composition": nucComp,
            "Nucleic Acid Count": sum(nucComp.values())
        }

        # Updates Total Amino Acid Composition
        for k, v in aaComp.items():
            self.genomeFeatures["Total Amino Acid Composition"][k] = \
            self.genomeFeatures["Total Amino Acid Composition"].get(k, 0) + v
        # Updates Total Codon Composition
        for k, v in codonComp.items():
            self.genomeFeatures["Total Codon Composition"][k] = \
            self.genomeFeatures["Total Codon Composition"].get(k, 0) + v
        # Updates Total Codon Count
        self.genomeFeatures["Total Codon Count"] += sum(codonComp.values())
        # Updates Total Nucleic Acid Composition
        for k, v in nucComp.items():
            self.genomeFeatures["Total Nucleic Acid Composition"][k] = \
            self.genomeFeatures["Total Nucleic Acid Composition"].get(k, 0) + v
        # Updates Total Nucleic Acid Count
        self.genomeFeatures["Total Nucleic Acid Count"] += sum(nucComp.values())
        # Updates root dictionary
        self.genomeFeatures["Individual Sequences"].update({})


    def getTotalAAComposition(self):
        """Retrieves the total amino acid composition across all analyzed
            sequences.

        Args:
            None

        Returns:
            The total amino acid composition.

        Raises:
            None
        """

        return self.genomeFeatures["Total Amino Acid Composition"]


    # Defines an alias for getTotalAAComposition
    aaComp = getTotalAAComposition

    def getTotalAACount(self, aa=None):
        """Retrieves the total count of a specific amino acid or all amino acids.

        Args:
            aa: The amino acid for which to retrieve the count.
                If None, the total count for all amino acids is returned.

        Returns:
            The total count of the specified amino acid or the sum of
                counts for all amino acids.

        Raises:
            None
        """

        if aa:
            return self.genomeFeatures["Total Amino Acid Composition"][aa]
        else:
            return sum(self.genomeFeatures["Total Amino Acid Composition"].values())


    # Defines an alias for getTotalAACount
    aaCount = getTotalAACount

    def getTotalAARelativeContent(self, aas):
        """Calculates the relative content of specified amino acids as a
            percentage of total amino acids.

        Args:
            aas (list[str] or str): A list of amino acids or a single amino acid.

        Returns:
            The relative content of specified amino acids as a percentage
                of the total.

        Raises:
            None
        """

        if not isinstance(aas, list):
            aas = [aas]
        aaCount = sum(self.genomeFeatures["Total Amino Acid Composition"][aa.upper()] \
                      for aa in aas)
        totalAAs = sum(self.genomeFeatures["Total Amino Acid Composition"].values())
        return (aaCount / totalAAs) * 100 if totalAAs > 0 else 0


    def getTotalCodonComposition(self):
        """Retrieves the total codon composition across all analyzed sequences.

        Returns:
            The total codon composition.

        Raises:
            None
        """

        return self.genomeFeatures["Total Codon Composition"]


    # Defines an alias for getTotalCodonComposition
    codonComp = getTotalCodonComposition

    def getTotalCodonCount(self, codon=None):
        """Retrieves the total count of a specific codon or all codons.

        Args:
            codon (str, optional): The codon for which to retrieve the count. If None,
                                the total count for all codons is returned.

        Returns:
            int: The total count of the specified codon or the sum of counts for all codons.

        Raises:
            None
        """

        if codon:
            return self.genomeFeatures["Total Codon Composition"][codon]
        else:
            return sum(self.genomeFeatures["Total Codon Composition"].values())


    # Defines an alias for getTotalCodonCount
    codonCount = getTotalCodonCount

    def getTotalCodonRelativeContent(self, codons):
        """Calculates the relative content of specified codons as a percentage
            of total codons.

        Args:
            codons: A list of codons or a single codon.

        Returns:
            The relative content of specified codons as a percentage of
                the total.
        """

        # Ensure codons is a list for uniform processing
        if not isinstance(codons, list):
            codons = [codons]

        # Calculate the total count of specified codons
        codonCount = sum(self.genomeFeatures["Total Codon Composition"].get(codon.upper(), 0) \
                         for codon in codons)

        # Calculate the sum of all codon counts
        totalCodons = sum(self.genomeFeatures["Total Codon Composition"].values())

        # Return the relative content as a percentage
        return (codonCount / totalCodons) * 100 if totalCodons > 0 else 0


    def getTotalNucComposition(self):
        """Retrieves the total nucleic acid composition across all analyzed sequences.

        Returns:
            The total nucleic acid composition.

        Raises:
            None
        """

        return self.genomeFeatures["Total Nucleic Acid Composition"]


    # Defines an alias for getTotalNucComposition
    nucComp = getTotalNucComposition

    def getTotalNucCount(self, nuc=None):
        """Calculates the GC content as a percentage of total nucleic acids.

        Returns:
            float: The GC content percentage.

        Raises:
            None
        """
        if nuc:
            return self.genomeFeatures["Total Nucleic Acid Composition"][nuc]
        else:
            return sum(self.genomeFeatures["Total Nucleic Acid Composition"].values())


    # Defines an alias for getTotalNucCount
    nucCount = getTotalNucCount

    def getTotalGCContent(self):
        """Retrieves the total count of a specific nucleotide or all nucleotides.

        Args:
            nuc: The nucleotide for which to retrieve the count. If None,
                                the total count for all nucleotides is returned.

        Returns:
            The total count of the specified nucleotide or the sum of counts for all nucleotides.

        Raises:
            None
        """

        gcCount = self.genomeFeatures["Total Nucleic Acid Composition"].get('G', 0)
        + self.genomeFeatures["Total Nucleic Acid Composition"].get('C', 0)
        totalNuc = sum(self.genomeFeatures["Total Nucleic Acid Composition"].values())
        return (gcCount / totalNuc) * 100 if totalNuc > 0 else 0


    def getTotalNucRelativeContent(self, nucs):
        """Calculates the relative content of specified nucleotides as a percentage of total nucleotides.

        Args:
            nucs: A list of nucleotides or a single nucleotide.

        Returns:
            The relative content of specified nucleotides as a percentage of the total.

        Raises:
            None
        """
        # Ensure nucs is a list for uniform processing
        if not isinstance(nucs, list):
            nucs = [nucs]

        # Calculate the total count of specified nucleotides
        nucCount = sum(self.genomeFeatures["Total Nucleic Acid Composition"].get(nuc.upper(), 0) \
                       for nuc in nucs)

        # Calculate the sum of all nucleotide counts
        totalNucs = sum(self.genomeFeatures["Total Nucleic Acid Composition"].values())

        # Return the relative content as a percentage
        return (nucCount / totalNucs) * 100 if totalNucs > 0 else 0
