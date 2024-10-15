import numpy as np
from agentBased.utils import rotMat


def angleMove(position, velocity, dt):
    theta = np.random.randn() * np.pi / 40
    velocity = np.matmul(rotMat(theta), velocity)
    position += dt * velocity
    return position, velocity


def cartesianMove(position, velocity, dt):
    vectors = [np.array([1, 0]), np.array([-1, 0]), np.array([0, 1]), np.array([0, -1])]
    v = np.linalg.norm(velocity)
    chosen_vector = v*vectors[np.random.choice(len(vectors))]
    position += chosen_vector * dt
    return position, velocity



class AndersonChaplain:
    def __init__(self, velocity, dt):
        self.h = np.linalg.norm(velocity)


    def andersonMove(self, position, velocity, dt):
        h = np.linalg.norm(velocity)*dt