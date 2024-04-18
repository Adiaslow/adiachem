"""
Module for calculating Mean Squared Logarithmic Error (MSLE).
"""

import numpy as np

def calculate_msle(reference_values, test_values):
    """
    Calculates the Mean Squared Logarithmic Error (MSLE) between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.

    Returns:
        float: The calculated MSLE value.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    log_squared_differences = np.square(np.log1p(reference_values) - np.log1p(test_values))
    msle = np.mean(log_squared_differences)

    return msle
