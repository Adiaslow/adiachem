"""
Module for calculating Explained Variance Score.
"""

import numpy as np

def calculate_explained_variance(reference_values, test_values):
    """
    Calculates the Explained Variance Score between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.

    Returns:
        float: The calculated Explained Variance Score.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    reference_variance = np.var(reference_values)
    residual_variance = np.var(reference_values - test_values)
    explained_variance = 1 - (residual_variance / reference_variance)

    return explained_variance
