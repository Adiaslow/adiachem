"""
Module for the AlignmentBase class, which provides a base class for sequence alignment algorithms.
"""

import numpy as np
import matplotlib.pyplot as plt

class AlignmentBase:
    """Base class for sequence alignment algorithms.

    Attributes:
        match_score (int): Score for matching elements.
        mismatch_score (int): Score for mismatching elements.
        gap_penalty (int): Score for gaps.
        score_matrix (numpy.ndarray): Matrix of scores for all possible alignments.

    Methods:
        align: Perform sequence alignment on two sequences.
        _traceback: Trace back through the score matrix to find the optimal alignment.
        _traceback_path: Trace back through the score matrix to find the optimal alignment path.
        visualize_alignment: Visualize the alignment of two sequences.
        visualize_alignment_matrix: Visualize the alignment matrix.
    """

    def _traceback_path(self, seq1, seq2, score_matrix):
        """Perform traceback to find the optimal alignment path.

        Args:
            seq1 (str): First sequence.
            seq2 (str): Second sequence.
            score_matrix (numpy.ndarray): Matrix of alignment scores.

        Returns:
            list: List of tuples representing the optimal alignment path.

        Raises:
            NotImplementedError: This method must be implemented by subclasses.
        """

        raise NotImplementedError("Subclasses must implement _traceback_path method")

    def visualize_alignment(self, align1, align2, width=60):
        """Visualize the alignment of two sequences.

        Args:
            align1 (str): First aligned sequence.
            align2 (str): Second aligned sequence.
            width (int): Width of each line in the visualization.

        Returns:
            None

        Raises:
            None
        """

        alignment_strings = []
        match_string = ""
        for a1, a2 in zip(align1, align2):
            if a1 == a2:
                match_string += "|"
            elif a1 == "-" or a2 == "-":
                match_string += " "
            else:
                match_string += "."

        num_lines = (len(align1) + width - 1) // width
        for i in range(num_lines):
            start = i * width
            end = min(start + width, len(align1))
            alignment_strings.append(f"Seq1: {align1[start:end]}")
            alignment_strings.append(f"      {match_string[start:end]}")
            alignment_strings.append(f"Seq2: {align2[start:end]}")
            alignment_strings.append("")

        alignment_string = "\n".join(alignment_strings)
        print(alignment_string)

    def visualize_alignment_matrix(self, seq1, seq2, score_matrix):
        """Visualize the alignment matrix.

        Args:
            seq1 (str): First sequence.
            seq2 (str): Second sequence.
            score_matrix (numpy.ndarray): Matrix of alignment scores.

        Returns:
            None

        Raises:
            None
        """

        fig, ax = plt.subplots(figsize=(10, 10))
        im = ax.imshow(score_matrix, cmap='viridis', origin='lower')

        ax.set_xticks(np.arange(len(seq2) + 1))
        ax.set_yticks(np.arange(len(seq1) + 1))
        ax.set_xticklabels([''] + list(seq2), fontsize=6)
        ax.set_yticklabels([''] + list(seq1), fontsize=6)

        ax.set_xlabel('Sequence 2', fontsize=10)
        ax.set_ylabel('Sequence 1', fontsize=10)
        ax.set_title('Alignment Matrix', fontsize=16)

        cbar = ax.figure.colorbar(im, ax=ax, shrink=0.77)
        cbar.ax.set_ylabel('Score', rotation=-90, va='bottom', fontsize=14)

        path = self._traceback_path(seq1, seq2, score_matrix)
        path_coords = [(p[1], p[0]) for p in path]
        path_coords.append((path_coords[-1][0] + 1, path_coords[-1][1] + 1))
        ax.plot(*zip(*path_coords[:-1]), color='red', linewidth=2, alpha=0.7)

        fig.tight_layout()
        plt.show()
