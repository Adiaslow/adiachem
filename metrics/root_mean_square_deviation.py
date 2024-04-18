"""
Module for calculating Root Mean Square Deviation (RMSD).
"""

import numpy as np

def calculate_rmsd(reference_coordinates, test_coordinates):
    """
    Calculates the Root Mean Square Deviation (RMSD) between two sets of coordinates.

    Args:
        reference_coordinates (numpy.ndarray): Reference coordinates.
        test_coordinates (numpy.ndarray): Test coordinates.

    Returns:
        float: The calculated RMSD value.
    """
    if reference_coordinates.shape != test_coordinates.shape:
        raise ValueError("Reference and test coordinates must have the same shape.")

    squared_differences = np.square(reference_coordinates - test_coordinates)
    sum_of_squares = np.sum(squared_differences, axis=1)
    mean_of_squares = np.mean(sum_of_squares)
    rmsd = np.sqrt(mean_of_squares)

    return rmsd
