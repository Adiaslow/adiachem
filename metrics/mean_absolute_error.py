"""
Module for calculating Mean Absolute Error (MAE).
"""

import numpy as np

def calculate_mae(reference_values, test_values):
    """
    Calculates the Mean Absolute Error (MAE) between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.

    Returns:
        float: The calculated MAE value.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    absolute_differences = np.abs(reference_values - test_values)
    mae = np.mean(absolute_differences)

    return mae
