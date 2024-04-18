"""
Module for calculating Mean Squared Error (MSE).
"""

import numpy as np

def calculate_mse(reference_values, test_values):
    """
    Calculates the Mean Squared Error (MSE) between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.

    Returns:
        float: The calculated MSE value.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    squared_differences = np.square(reference_values - test_values)
    mse = np.mean(squared_differences)

    return mse
