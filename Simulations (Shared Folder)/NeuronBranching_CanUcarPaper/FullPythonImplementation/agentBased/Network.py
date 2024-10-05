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

        for ducts in self.ductsToBeAdded: ## The newly created ducts (due to branching) are not added to the
            # duct list right away as it can lead to buggy behaviour (since the loop is happening on the duct list and
            # it should not be altered)
            self.ductsList.append(ducts)


        self.ductsToBeAdded = []  ## To make sure that the ducts are added and the list is empty for the next round

    def handleBranching(self, initPos):
        v = np.array([1.0,0.0])
        child1_angle = np.random.randn()*np.pi/100.0 + np.pi/8
        child2_angle = np.random.randn()*np.pi/100.0 - np.pi/8
        self.ductsToBeAdded.append(Duct(initPos,np.matmul(rotMat(child1_angle),v),self))
        self.ductsToBeAdded.append(Duct(initPos,np.matmul(rotMat(child2_angle),v),self))


