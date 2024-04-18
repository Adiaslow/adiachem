"""
Module for calculating Root Mean Square Error (RMSE).
"""

import numpy as np

def calculate_rmse(reference_values, test_values):
    """
    Calculates the Root Mean Square Error (RMSE) between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.

    Returns:
        float: The calculated RMSE value.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    squared_differences = np.square(reference_values - test_values)
    mean_of_squares = np.mean(squared_differences)
    rmse = np.sqrt(mean_of_squares)

    return rmse
