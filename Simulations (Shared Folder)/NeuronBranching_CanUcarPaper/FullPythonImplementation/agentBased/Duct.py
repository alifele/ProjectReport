from Tip import Tip
import numpy as np

class Duct:
    def __init__(self,initPos, initVel, network):
        self.tip = Tip(initPos,initVel, self)
        self.tailCurve = []
        self.tailCurve.append(self.tip.position.tolist())
        self.network = network

    def update(self):
        if self.tip.isAlive:
            self.tip.update()
            self.tailCurve.append(self.tip.position.tolist())

