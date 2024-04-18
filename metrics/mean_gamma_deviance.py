"""
Module for calculating Mean Gamma Deviance.
"""

import numpy as np

def calculate_mean_gamma_deviance(reference_values, test_values):
    """
    Calculates the Mean Gamma Deviance between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.

    Returns:
        float: The calculated Mean Gamma Deviance.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    gamma_deviance = 2 * (np.log(reference_values / test_values) + (test_values - reference_values) / reference_values)
    mean_gamma_deviance = np.mean(gamma_deviance)

    return mean_gamma_deviance
