"""
Module for calculating Pinball Loss.
"""

import numpy as np

def calculate_pinball_loss(reference_values, test_values, quantile):
    """
    Calculates the Pinball Loss between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.
        quantile (float): Quantile parameter for the Pinball Loss.

    Returns:
        float: The calculated Pinball Loss.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    residuals = reference_values - test_values
    pinball_loss = np.mean(np.where(residuals >= 0, quantile * residuals, (quantile - 1) * residuals))

    return pinball_loss
