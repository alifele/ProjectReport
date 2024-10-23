import numpy as np

from agentBased.Data import Data
from agentBased.utils import rotMat


def angleMove(position, velocity, dt):

    ## Influenced by external vector field
    thetaDifference = angle_between(vectorField(position),velocity)


    # DeltaTheta = np.random.randn() * np.pi / 40
    DeltaTheta = 2.5*np.random.randn()*np.pi*dt - 2*thetaDifference*dt
    velocity = np.matmul(rotMat(DeltaTheta), velocity)

    position += dt * velocity
    return position, velocity



def cartesianMove(position, velocity, dt):
    vectors = [np.array([1, 0]), np.array([-1, 0]), np.array([0, 1]), np.array([0, -1])]
    v = np.linalg.norm(velocity)
    chosen_vector = v*vectors[np.random.choice(len(vectors))]
    position += chosen_vector * dt
    return position, velocity


def vectorField(position):
    # abs = np.linalg.norm(position)
    x = position[0]
    y = position[1]

    ### Spiral
    # vx = 0.1 -y + 0.4*x
    # vy = 0.1 + x + 0.4*y

    ### Simple
    # vx = 1.0
    # vy = 1.0

    ### Funny
    # vx = x - x*y + 0.1
    # vy = x-y + 0.1


    ### Tumor
    x, y = position
    x0, y0 = Data.tumorPos
    R = Data.tumorRad
    exp_term = np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2 * R ** 2))
    vx = (x - x0) / R ** 2 * exp_term
    vy = (y - y0) / R ** 2 * exp_term

    # Return the gradient as a numpy array
    return -np.array([vx, vy])


class AndersonChaplain:
    def __init__(self, velocity, dt):
        self.h = np.linalg.norm(velocity)


    def andersonMove(self, position, velocity, dt):
        h = np.linalg.norm(velocity)*dt


def angle_between(v1, v2):
    # Ensure the vectors are numpy arrays
    v1 = np.array(v1)
    v2 = np.array(v2)

    # Normalize the vectors
    v1_normalized = v1 / np.linalg.norm(v1)
    v2_normalized = v2 / np.linalg.norm(v2)

    # Calculate the dot product
    dot_product = np.dot(v1_normalized, v2_normalized)

    # Calculate the determinant (for the direction of the angle)
    determinant = np.linalg.det([v1_normalized, v2_normalized])

    # Compute the angle using arctan2 for range [-pi, pi]
    angle = np.arctan2(determinant, dot_product)

    return angle
