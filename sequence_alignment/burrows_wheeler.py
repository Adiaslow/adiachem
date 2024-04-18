"""
Module for performing BWA (Burrows-Wheeler Aligner) on sequences.
"""

from alignment_base import AlignmentBase
import numpy as np

class BWA(AlignmentBase):
    """Class for performing BWA on sequences.

    Attributes:
        reference (str): Reference sequence.
        bwt (str): Burrows-Wheeler Transform of the reference sequence.
        suffix_array (numpy.ndarray): Suffix array of the reference sequence.

    Methods:
        _build_bwt: Build the Burrows-Wheeler Transform of the sequence.
        _build_suffix_array: Build the suffix array of the sequence.
        _backward_search: Perform backward search to find the range of suffixes that match the pattern.
        align: Align the query sequence against the reference sequence.
    """

    def __init__(self, reference):
        """Initialize the BWA object.

        Args:
            reference (str): Reference sequence.

        Returns:
            None

        Raises:
            None
        """
        self.reference = reference
        self.bwt = self._build_bwt(reference)
        self.suffix_array = self._build_suffix_array(reference)

    def _build_bwt(self, sequence):
        """Build the Burrows-Wheeler Transform of the sequence.

        Args:
            sequence (str): Input sequence.

        Returns:
            str: Burrows-Wheeler Transform of the sequence.

        Raises:
            None
        """
        sequence += '$'
        rotations = [sequence[i:] + sequence[:i] for i in range(len(sequence))]
        rotations.sort()
        bwt = ''.join([r[-1] for r in rotations])
        return bwt

    def _build_suffix_array(self, sequence):
        """Build the suffix array of the sequence.

        Args:
            sequence (str): Input sequence.

        Returns:
            numpy.ndarray: Suffix array of the sequence.

        Raises:
            None
        """
        sequence += '$'
        suffixes = [(sequence[i:], i) for i in range(len(sequence))]
        suffixes.sort()
        suffix_array = np.array([i for _, i in suffixes])
        return suffix_array

    def _backward_search(self, pattern):
        """Perform backward search to find the range of suffixes that match the pattern.

        Args:
            pattern (str): Pattern to search for.

        Returns:
            tuple: Range of suffixes that match the pattern (start, end).

        Raises:
            None
        """
        start, end = 0, len(self.bwt) - 1
        for c in reversed(pattern):
            start = self.bwt[:start].count(c) + 1
            end = self.bwt[:end + 1].count(c)
            if start > end:
                return None
        return start, end

    def align(self, query):
        """Align the query sequence against the reference sequence.

        Args:
            query (str): Query sequence.

        Returns:
            list: List of tuples representing the alignments (position, cigar).

        Raises:
            None
        """
        alignments = []
        for i in range(len(query)):
            pattern = query[i:]
            match = self._backward_search(pattern)
            if match:
                start, end = match
                positions = self.suffix_array[start:end + 1]
                for pos in positions:
                    cigar = self._compute_cigar(query, self.reference[pos:])
                    alignments.append((pos, cigar))
        return alignments

    def _compute_cigar(self, query, reference):
        """Compute the CIGAR string for the alignment.

        Args:
            query (str): Query sequence.
            reference (str): Reference sequence.

        Returns:
            str: CIGAR string representing the alignment.

        Raises:
            None
        """
        cigar = ''
        i, j = 0, 0
        while i < len(query) and j < len(reference):
            if query[i] == reference[j]:
                cigar += 'M'
                i += 1
                j += 1
            elif query[i] != reference[j]:
                cigar += 'X'
                i += 1
                j += 1
            elif i == len(query):
                cigar += 'D'
                j += 1
            elif j == len(reference):
                cigar += 'I'
                i += 1
        if i < len(query):
            cigar += 'I' * (len(query) - i)
        elif j < len(reference):
            cigar += 'D' * (len(reference) - j)
        return cigar
