import numpy as np
from dask.array import angle

from agentBased.MovementModels import *


class Tip:
    def __init__(self, initPos, initVel, duct):
        self.position = np.copy(initPos)
        self.velocity = np.copy(initVel) ## The absolute value of the velocity
        self.collisionData = None ## this will be replaced with collision information if collides: [Stopperduct, stopperPoint]
        self.isAlive = True
        self.duct = duct
        self.dt = self.duct.network.dt

    def update(self):
        if self.collisionData == None: ## Is not colliding
            self.move()
            self.branch()
        else:
            stopperDuct = self.collisionData[0]
            stopperPoint = self.collisionData[1]
            stopperDuct.collisionList.append(self)
            self.position = stopperPoint
            self.isAlive = False


    def move(self):
        # self.position, self.velocity = cartesianMove(self.position, self.velocity, self.dt)
        self.position, self.velocity = angleMove(self.position, self.velocity, self.dt)


    def branch(self):
        if np.random.random() < 0.022:
            self.duct.network.handleBranching(self.position, self.duct)
            self.isAlive = False
