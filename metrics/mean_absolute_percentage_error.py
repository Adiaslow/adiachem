"""
Module for calculating Mean Absolute Percentage Error (MAPE).
"""

import numpy as np

def calculate_mape(reference_values, test_values):
    """
    Calculates the Mean Absolute Percentage Error (MAPE) between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.

    Returns:
        float: The calculated MAPE value.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    absolute_percentage_errors = np.abs((reference_values - test_values) / reference_values)
    mape = np.mean(absolute_percentage_errors) * 100

    return mape
