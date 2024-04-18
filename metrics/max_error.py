"""
Module for calculating Max Error.
"""

import numpy as np

def calculate_max_error(reference_values, test_values):
    """
    Calculates the Max Error between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.

    Returns:
        float: The calculated Max Error value.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    absolute_errors = np.abs(reference_values - test_values)
    max_error = np.max(absolute_errors)

    return max_error
