"""
Module for calculating Median Absolute Error (MedAE).
"""

import numpy as np

def calculate_medae(reference_values, test_values):
    """
    Calculates the Median Absolute Error (MedAE) between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.

    Returns:
        float: The calculated MedAE value.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    absolute_errors = np.abs(reference_values - test_values)
    medae = np.median(absolute_errors)

    return medae
