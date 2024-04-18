"""
Module for finding the optimal translation vector and rotation matrix to superimpose two arrays of coordinates using Singular Value Decomposition (SVD).
"""
from numpy import array, mean
from numpy.linalg import svd

def singular_value_decomposition(coords1, coords2):
    """Finds optimal translation vector and rotation matrix to superimpose two arrays of coordinates using Singular Value Decomposition (SVD).

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

    _, _, vt = svd(centered_coords2.T @ centered_coords1)

    rotation_matrix = vt.T
    translation_vector = centroid1 - centroid2 @ rotation_matrix

    return translation_vector, rotation_matrix
