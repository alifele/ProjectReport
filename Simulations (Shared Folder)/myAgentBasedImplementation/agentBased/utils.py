import numpy as np

def rotMat(theta):
    """
    Returns the 2D rotation matrix for a given angle theta (in radians).

    Parameters:
    theta (float): The rotation angle in radians.

    Returns:
    numpy.ndarray: A 2x2 rotation matrix.
    """
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])