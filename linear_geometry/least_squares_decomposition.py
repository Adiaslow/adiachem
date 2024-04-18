"""
Module for finding the optimal translation vector and rotation matrix to superimpose two arrays of coordinates using Least Squares Fitting (LSF).
"""
from numpy import array, mean
from scipy.linalg import orthogonal_procrustes

def least_squares_decomposition(coords1, coords2):
    """Finds optimal translation vector and rotation matrix to superimpose two arrays of coordinates using Least Squares Fitting (LSF).

    Args:
        coords1 (np.array): Coordinates of the first set.
        coords2 (np.array): Coordinates of the second set.

    Returns:
        tuple: Translation vector and rotation matrix.

    Raises:
        None
    """
    coords1 = array(coords1)
    coords2 = array(coords2)

    centroid1 = mean(coords1, axis=0)
    centroid2 = mean(coords2, axis=0)

    centered_coords1 = coords1 - centroid1
    centered_coords2 = coords2 - centroid2

    rotation_matrix, _ = orthogonal_procrustes(centered_coords2, centered_coords1)
    translation_vector = centroid1 - centroid2 @ rotation_matrix

    return translation_vector, rotation_matrix
