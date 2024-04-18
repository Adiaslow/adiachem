"""
Module for performing Smith-Waterman local alignment on sequences.
"""

from alignment_base import AlignmentBase
import numpy as np

class SmithWaterman(AlignmentBase):
    """Class for performing Smith-Waterman local alignment on sequences.

    Attributes:
        match_score (int): Score for matching characters.
        mismatch_score (int): Score for mismatching characters.
        gap_penalty (int): Penalty for gaps.

    Methods:
        align: Perform Smith-Waterman local alignment on the input sequences.
        _traceback: Perform traceback to obtain the aligned sequences.
        _traceback_path: Perform traceback to obtain the aligned sequences and the alignment path.
    """

    def __init__(self, match_score=2, mismatch_score=-1, gap_penalty=-1):
        """Initialize the SmithWaterman object.

        Args:
            match_score (int): Score for matching characters.
            mismatch_score (int): Score for mismatching characters.
            gap_penalty (int): Penalty for gaps.

        Returns:
            None

        Raises:
            None
        """
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty

    def align(self, seq1, seq2):
        """Perform Smith-Waterman local alignment on the input sequences.

        Args:
            seq1 (str): First input sequence.
            seq2 (str): Second input sequence.

        Returns:
            tuple: Tuple containing the aligned sequences and the alignment score.
        """
        m, n = len(seq1), len(seq2)
        matrix = np.zeros((m+1, n+1), dtype=int)

        max_score = 0
        max_pos = None

        for i in range(1, m+1):
            for j in range(1, n+1):
                score_diag = matrix[i-1][j-1] + (self.match_score if seq1[i-1] == seq2[j-1] else self.mismatch_score)
                score_up = matrix[i-1][j] + self.gap_penalty
                score_left = matrix[i][j-1] + self.gap_penalty
                matrix[i][j] = max(0, score_diag, score_up, score_left)
                if matrix[i][j] > max_score:
                    max_score = matrix[i][j]
                    max_pos = (i, j)

        aligned_seq1, aligned_seq2 = self._traceback(seq1, seq2, matrix, max_pos)
        return aligned_seq1, aligned_seq2, max_score

    def _traceback(self, seq1, seq2, matrix, pos):
        """Perform traceback to obtain the aligned sequences.

        Args:
            seq1 (str): First input sequence.
            seq2 (str): Second input sequence.
            matrix (numpy.ndarray): Scoring matrix.
            pos (tuple): Position of the maximum score in the matrix.

        Returns:
            tuple: Tuple containing the aligned sequences.
        """
        aligned_seq1, aligned_seq2 = [], []
        i, j = pos
        while matrix[i][j] != 0:
            if matrix[i][j] == matrix[i-1][j-1] + (self.match_score if seq1[i-1] == seq2[j-1] else self.mismatch_score):
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif matrix[i][j] == matrix[i-1][j] + self.gap_penalty:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1

        aligned_seq1 = ''.join(reversed(aligned_seq1))
        aligned_seq2 = ''.join(reversed(aligned_seq2))
        return aligned_seq1, aligned_seq2

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
