import numpy as np
from utils import rotMat


class Tip:
    def __init__(self, initPos, initVel, duct):
        self.position = np.copy(initPos)
        self.velocity = np.copy(initVel) ## The absolute value of the velocity
        self.isAlive = True
        self.duct = duct

        self.dt = 0.01

    def update(self):
        self.move()
        self.branch()



    def move(self):
        theta = np.random.randn() * np.pi / 50
        self.velocity = np.matmul(rotMat(theta), self.velocity)
        self.position += self.dt * self.velocity

    def branch(self):
        if np.random.random() < 0.01:
            self.duct.network.handleBranching(np.copy(self.position))
            self.isAlive = False

