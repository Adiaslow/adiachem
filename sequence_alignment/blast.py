"""
Module for performing BLAST (Basic Local Alignment Search Tool) on sequences.
"""

from alignment_base import AlignmentBase
from collections import defaultdict
from itertools import product

class BLAST(AlignmentBase):
    """Class for performing BLAST on sequences.

    Attributes:
        word_size (int): Size of the word for BLAST algorithm.
        score_matrix (dict): Scoring matrix for BLAST alignment.
        threshold (int): Threshold score for BLAST hits.

    Methods:
        _generate_score_matrix: Generate the score matrix for BLAST alignment.
        _find_hits: Find high-scoring segment pairs (HSPs) between the query and database sequences.
        _generate_words: Generate all possible words of the specified size from the sequence.

    Raises:
        None
    """

    def __init__(self, word_size=3, match_score=1, mismatch_score=-1, threshold=11):
        """Initialize the BLAST object.

        Args:
            word_size (int): Size of the word for BLAST algorithm.
            match_score (int): Score for matching elements.
            mismatch_score (int): Score for mismatching elements.
            threshold (int): Threshold score for BLAST hits.

        Returns:
            None

        Raises:
            None
        """
        self.word_size = word_size
        self.score_matrix = self._generate_score_matrix(match_score, mismatch_score)
        self.threshold = threshold

    def _generate_score_matrix(self, match_score, mismatch_score):
        """Generate the score matrix for BLAST alignment.

        Args:
            match_score (int): Score for matching elements.
            mismatch_score (int): Score for mismatching elements.

        Returns:
            dict: Scoring matrix for BLAST alignment.

        Raises:
            None
        """
        score_matrix = defaultdict(lambda: defaultdict(int))
        for c1, c2 in product('ACDEFGHIKLMNPQRSTVWY', repeat=2):
            if c1 == c2:
                score_matrix[c1][c2] = match_score
            else:
                score_matrix[c1][c2] = mismatch_score
        return score_matrix

    def _find_hits(self, query, database):
        """Find high-scoring segment pairs (HSPs) between the query and database sequences.

        Args:
            query (str): Query sequence.
            database (str): Database sequence.

        Returns:
            list: List of tuples representing the HSPs.

        Raises:
            None
        """
        hits = []
        query_words = self._generate_words(query)
        database_words = self._generate_words(database)
        for i, word in enumerate(query_words):
            if word in database_words:
                for j in database_words[word]:
                    score = self._extend_hit(query, database, i, j)
                    if score >= self.threshold:
                        hits.append((i, j, score))
        return hits

    def _generate_words(self, sequence):
        """Generate all possible words of the specified size from the sequence.

        Args:
            sequence (str): Input sequence.

        Returns:
            dict: Dictionary of words and their positions in the sequence.

        Raises:
            None
        """
        words = defaultdict(list)
        for i in range(len(sequence) - self.word_size + 1):
            word = sequence[i:i + self.word_size]
            words[word].append(i)
        return words

    def _extend_hit(self, query, database, query_pos, database_pos):
        """Extend the hit in both directions to find the highest-scoring segment pair.

        Args:
            query (str): Query sequence.
            database (str): Database sequence.
            query_pos (int): Starting position of the hit in the query sequence.
            database_pos (int): Starting position of the hit in the database sequence.

        Returns:
            int: Score of the highest-scoring segment pair.

        Raises:
            None
        """
        score = 0
        query_left, query_right = query_pos, query_pos + self.word_size
        database_left, database_right = database_pos, database_pos + self.word_size

        while query_left > 0 and database_left > 0 and \
                self.score_matrix[query[query_left - 1]][database[database_left - 1]] > 0:
            score += self.score_matrix[query[query_left - 1]][database[database_left - 1]]
            query_left -= 1
            database_left -= 1

        while query_right < len(query) and database_right < len(database) and \
                self.score_matrix[query[query_right]][database[database_right]] > 0:
            score += self.score_matrix[query[query_right]][database[database_right]]
            query_right += 1
            database_right += 1

        return score

    def search(self, query, database):
        """Perform BLAST search on the query sequence against the database sequence.

        Args:
            query (str): Query sequence.
            database (str): Database sequence.

        Returns:
            list: List of tuples representing the high-scoring segment pairs (HSPs).

        Raises:
            None
        """
        hits = self._find_hits(query, database)
        return hits

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
