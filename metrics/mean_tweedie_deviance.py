"""
Module for calculating Mean Tweedie Deviance.
"""

import numpy as np

def calculate_mean_tweedie_deviance(reference_values, test_values, power):
    """
    Calculates the Mean Tweedie Deviance between two sets of values.

    Args:
        reference_values (numpy.ndarray): Reference values.
        test_values (numpy.ndarray): Test values.
        power (float): Power parameter for the Tweedie distribution.

    Returns:
        float: The calculated Mean Tweedie Deviance.
    """

    if reference_values.shape != test_values.shape:
        raise ValueError("Reference and test values must have the same shape.")

    if power == 0:
        tweedie_deviance = reference_values * np.log(reference_values / test_values) - reference_values + test_values
    elif power == 1:
        tweedie_deviance = reference_values * np.log(reference_values / test_values) - reference_values + test_values
    else:
        tweedie_deviance = (reference_values ** (2 - power) - (1 - power) * reference_values * test_values ** (1 - power)) / (1 - power) - (reference_values ** (2 - power) - (1 - power) * reference_values ** 2) / (1 - power)

    mean_tweedie_deviance = np.mean(tweedie_deviance)

    return mean_tweedie_deviance
