"""
Module for calculating Mean Poisson Deviance.
"""

import numpy as np

def calculate_mean_poisson_deviance(reference_values, test_values):
    """
    Calculates the Mean Poisson Deviance between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.

    Returns:
        float: The calculated Mean Poisson Deviance.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    poisson_deviance = 2 * (reference_values * np.log(reference_values / test_values) - reference_values + test_values)
    mean_poisson_deviance = np.mean(poisson_deviance)

    return mean_poisson_deviance
