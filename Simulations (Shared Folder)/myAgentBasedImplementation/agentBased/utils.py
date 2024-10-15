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


def calculate_curve_length(points):
    """
    Calculate the length of a curve given its points as a list or array of coordinates.

    Parameters:
    points (list or array): List of points (each point can be in 2D or 3D).

    Returns:
    float: The length of the curve.
    """

    points = np.array(points)  # Ensure the points are a NumPy array for vectorized operations

    if points.shape[0] < 2:
        return 0.01

    # Compute the pairwise differences between consecutive points
    diffs = np.diff(points, axis=0)
    # Compute the Euclidean distances between consecutive points
    distances = np.linalg.norm(diffs, axis=1)
    # Sum up all the distances to get the total length
    return np.sum(distances)
