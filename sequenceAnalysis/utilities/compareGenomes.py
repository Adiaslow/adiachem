class CompareGenomes:
    """A class for comparing two genomes.

    Attributes:
        genome1: A dictionary containing the first genome.
        genome2: A dictionary containing the second genome.

    Methods:
        compareGCContent: Compares the GC content of two genomes.
        compareAAComposition: Compares the amino acid composition of two
            genomes.
        compareRelativeCodonBias: Compares the relative codon bias of two
            genomes based on their total codon compositions.
    """

    def __init__(self, genome1, genome2):
        """The constructor for CompareGenomes.

        Args:
            genome1: A dictionary containing the first genome.
            genome2: A dictionary containing the second genome.

        Returns:
            None

        Raises:
            None
        """

        self.genome1 = genome1
        self.genome2 = genome2

    def compareGCContent(self):
        """Compares the GC content of two genomes.

        Args:
            None

        Returns:
            None

        Raises:
            None
        """

        gcContent1 = self.genome1.getTotalNucRelativeContent(['G', 'C'])
        gcContent2 = self.genome2.getTotalNucRelativeContent(['G', 'C'])
        print("\n          GC Content Comparison")
        print("-" * 40)
        print(f"           | Genome 1 | Genome 2 | Ratio")
        print("-" * 40)
        print(f"GC content |   {gcContent1:6.2f} "
            + f"|   {gcContent2:6.2f} "
            + f"|  {gcContent1 / gcContent2:.2f}")


    def compareAAComposition(self):
        """Compares the amino acid composition of two genomes.

        Args:
            None

        Returns:
            None

        Raises:
            None
        """

        aaComposition1 = self.genome1.getTotalAAComposition()
        aaComposition2 = self.genome2.getTotalAAComposition()
        print("\n   Amino Acid Composition Comparison")
        print("-" * 40)
        print("AA         | Genome 1 | Genome 2 | Ratio")
        print("-" * 40)
        for aa in sorted(aaComposition1.keys()):
            print(f"{aa}          |   {aaComposition1[aa]:6d} "
            + f"|   {aaComposition2[aa]:6d} "
            + f"|  {aaComposition1[aa]/aaComposition2[aa]:.2f}")


    def compareRelativeCodonBias(self):
        """Compares the relative codon bias of two genomes based on their total codon compositions.

        Args:
            None

        Returns:
            None

        Raises:
            None
        """
        print("\n     Relative Codon Bias Comparison")
        print("-" * 40)
        print("Codon      | Genome 1 | Genome 2 | Ratio")
        print("-" * 40)

        codonComp1 = self.genome1.genomeFeatures["Total Codon Composition"]
        codonComp2 = self.genome2.genomeFeatures["Total Codon Composition"]
        allCodons = set(codonComp1.keys()) | set(codonComp2.keys())
        totalCodons1 = sum(codonComp1.values())
        totalCodons2 = sum(codonComp2.values())

        for codon in sorted(allCodons):
            # Compute relative content for each genome
            relContent1 = (codonComp1.get(codon, 0) / totalCodons1) \
                * 100 if totalCodons1 > 0 else 0
            relContent2 = (codonComp2.get(codon, 0) / totalCodons2) \
                * 100 if totalCodons2 > 0 else 0

            # Compute ratio of relative contents, avoiding division by zero
            ratio = (relContent1 / relContent2) if relContent2 > 0 else 'N/A'
            print(f"{codon}        | {relContent1:8.2f}"
            + f" | {relContent2:8.2f} | "
            + (f" {ratio:3.2f}" if isinstance(ratio, float) else "N/A"))
