"""
Module for calculating Coefficient of Determination (R²).
"""

import numpy as np

def calculate_r2(reference_values, test_values):
    """
    Calculates the Coefficient of Determination (R²) between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.

    Returns:
        float: The calculated R² value.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    reference_mean = np.mean(reference_values)
    total_sum_of_squares = np.sum(np.square(reference_values - reference_mean))
    residual_sum_of_squares = np.sum(np.square(reference_values - test_values))
    r2 = 1 - (residual_sum_of_squares / total_sum_of_squares)

    return r2
