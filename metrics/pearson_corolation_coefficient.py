"""
Module for calculating Pearson Correlation Coefficient (PCC).
"""

import numpy as np

def calculate_pcc(reference_values, test_values):
    """
    Calculates the Pearson Correlation Coefficient (PCC) between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.

    Returns:
        float: The calculated PCC value.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    reference_mean = np.mean(reference_values)
    test_mean = np.mean(test_values)

    reference_std = np.std(reference_values)
    test_std = np.std(test_values)

    covariance = np.mean((reference_values - reference_mean) * (test_values - test_mean))
    pcc = covariance / (reference_std * test_std)

    return pcc
