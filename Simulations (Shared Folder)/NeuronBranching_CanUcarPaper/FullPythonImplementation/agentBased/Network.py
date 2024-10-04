from Duct import Duct
import numpy as np
from utils import rotMat

class Network:
    def __init__(self):
        self.ductsList = []
        self.ductsList.append(Duct(np.array([0.0,0.0]),np.array([1.0,0.0]), self))

        self.ductsToBeAdded = []


    def update(self):
        for ducts in self.ductsList:
            ducts.update()

        for ducts in self.ductsToBeAdded:
            self.ductsList.append(ducts)

        self.ductsToBeAdded = []

    def handleBranching(self, initPos):
        v = np.array([1.0,0.0])
        child1_angle = np.random.uniform(np.pi/4, np.pi/2)
        child2_angle = np.random.uniform(-np.pi/2, -np.pi/4)

        self.ductsToBeAdded.append(Duct(initPos,np.matmul(rotMat(child1_angle),v),self))
        self.ductsToBeAdded.append(Duct(initPos,np.matmul(rotMat(child2_angle),v),self))
        # self.ductsToBeAdded.append(Duct(np.array([0.0,0.0]),np.matmul(rotMat(child1_angle),v),self))
        # self.ductsToBeAdded.append(Duct(np.array([0.0,0.0]),np.matmul(rotMat(child2_angle),v),self))

