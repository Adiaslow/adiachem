"""
Module for aligning sequences using the Needleman-Wunsch algorithm.
"""

from alignment_base import AlignmentBase
import numpy as np
import matplotlib.pyplot as plt

class NeedlemanWunsch(AlignmentBase):
    """Class for aligning sequences using the Needleman-Wunsch algorithm.

    Attributes:
        match_score (int): Score for matching elements.
        mismatch_score (int): Score for mismatching elements.
        gap_penalty (int): Score for gaps.
        score_matrix (numpy.ndarray): Matrix of scores for all possible alignments.

    Methods:
        align: Perform Needleman-Wunsch alignment on two sequences.
        _traceback: Trace back through the score matrix to find the optimal alignment.
        _traceback_path: Trace back through the score matrix to find the optimal alignment path.

    """

    def __init__(self, match_score=1, mismatch_score=-1, gap_penalty=-1):
        """Initialize the NeedlemanWunsch object.

        Args:
            match_score (int): Score for matching elements.
            mismatch_score (int): Score for mismatching elements.
            gap_penalty (int): Score for gaps.

        Returns:
            None

        Raises:
            None
        """
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
        self.score_matrix = None

    def align(self, seq1, seq2):
        """Perform Needleman-Wunsch alignment on two sequences.

        Args:
            seq1 (str): First sequence to align.
            seq2 (str): Second sequence to align.

        Returns:
            tuple: Aligned sequence 1 and aligned sequence 2.
        """
        n, m = len(seq1), len(seq2)
        self.score_matrix = np.zeros((n + 1, m + 1))

        for i in range(1, n + 1):
            self.score_matrix[i, 0] = i * self.gap_penalty
        for j in range(1, m + 1):
            self.score_matrix[0, j] = j * self.gap_penalty

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                match = self.score_matrix[i - 1, j - 1] + (
                    self.match_score if seq1[i - 1] == seq2[j - 1] else self.mismatch_score
                )
                delete = self.score_matrix[i - 1, j] + self.gap_penalty
                insert = self.score_matrix[i, j - 1] + self.gap_penalty
                self.score_matrix[i, j] = max(match, delete, insert)

        align1, align2 = self._traceback(seq1, seq2, self.score_matrix)
        return align1, align2

    def _traceback(self, seq1, seq2, score_matrix):
        """Perform traceback to find the optimal alignment.

        Args:
            seq1 (str): First sequence to align.
            seq2 (str): Second sequence to align.
            score_matrix (numpy.ndarray): Matrix of scores for all possible alignments.

        Returns:
            tuple: Aligned sequence 1 and aligned sequence 2.
        """
        align1, align2 = "", ""
        i, j = len(seq1), len(seq2)
        while i > 0 and j > 0:
            if score_matrix[i, j] == score_matrix[i - 1, j - 1] + (
                self.match_score if seq1[i - 1] == seq2[j - 1] else self.mismatch_score
            ):
                align1 = seq1[i - 1] + align1
                align2 = seq2[j - 1] + align2
                i -= 1
                j -= 1
            elif score_matrix[i, j] == score_matrix[i - 1, j] + self.gap_penalty:
                align1 = seq1[i - 1] + align1
                align2 = "-" + align2
                i -= 1
            else:
                align1 = "-" + align1
                align2 = seq2[j - 1] + align2
                j -= 1

        while i > 0:
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        while j > 0:
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

        return align1, align2

    def _traceback_path(self, seq1, seq2, score_matrix):
        """Perform traceback to find the optimal alignment path.

        Args:
            seq1 (str): First sequence to align.
            seq2 (str): Second sequence to align.
            score_matrix (numpy.ndarray): Matrix of scores for all possible alignments.

        Returns:
            tuple: Aligned sequence 1 and aligned sequence 2.

        Raises:
            NotImplementedError: This method is not yet implemented.
        """
        raise NotImplementedError("This method is not yet implemented.")
