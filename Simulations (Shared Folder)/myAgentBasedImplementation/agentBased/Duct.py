from Tip import Tip
import numpy as np

class Duct:
    def __init__(self,initPos, initVel, network, parent):
        self.network = network
        self.parent = parent
        self.tailCurve = []
        self.collisionList = [] # Contains the tips that collided with this duct
        self.tip = Tip(initPos, initVel, self)
        self.tailCurve.append(self.tip.position.tolist())

    def update(self):
        if self.tip.isAlive:
            self.tip.update()
            self.tailCurve.append(self.tip.position.tolist())

