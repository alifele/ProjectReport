from numpy import number

from Duct import Duct
import numpy as np
from utils import rotMat
import matplotlib.pyplot as plt

from scipy.spatial import KDTree

class Network:
    def __init__(self):
        self.ductsList = []
        self.ductsToBeAdded = []  ## will add the newly generated ducts (due to branching) after when iteration on the
        # ductsList to update the ducts is completed.
        self.dt = 0.01
        self.vNorm = 1.0 # the absolute value of the velocity. Note that since the velocity will be changing by a rotation
        # matrix, then its norm will stay the same throughout
        self.pointToDuct = {}
        self.r_c = (self.vNorm * self.dt) * 1 # the r collision should be larger than the maximum step size of the walkers

        self.initialize()

    def initialize(self):
        self.ductsList.append(Duct(np.array([0.0,0.0]),np.array([self.vNorm,0.0]), self, None))
        self.ductsList.append(Duct(np.array([10.0,0.0]),np.array([0.0, self.vNorm]), self, None))
        self.ductsList.append(Duct(np.array([0.0,10.0]),np.array([0.0, -self.vNorm]), self, None))
        self.ductsList.append(Duct(np.array([10.0,10.0]),np.array([-self.vNorm, 0.0]), self, None))


    def update(self):
        self.collisionDetection()
        for ducts in self.ductsList:
            ducts.update()

        for ducts in self.ductsToBeAdded: ## The newly created ducts (due to branching) are not added to the
            # duct list right away as it can lead to buggy behaviour (since the loop is happening on the duct list and
            # it should not be altered)
            self.ductsList.append(ducts)



        self.ductsToBeAdded = []  ## To make sure that the ducts are added and the list is empty for the next round

    def handleBranching(self, initPos, parent):
        v = parent.tip.velocity
        # child1_angle = np.random.randn()*np.pi/100.0 + np.pi/8
        # child2_angle = np.random.randn()*np.pi/100.0 - np.pi/8
        child1_angle = np.pi/8
        child2_angle = -np.pi/8
        self.ductsToBeAdded.append(Duct(initPos,np.matmul(rotMat(child1_angle),v),self, parent))
        self.ductsToBeAdded.append(Duct(initPos,np.matmul(rotMat(child2_angle),v),self, parent))


    def collisionDetection(self):
        """
        ** This function will detect the collisions that will happen if we take the next steps for tips.
        For this we are using the fact that the radius of the collision detection should be larger than the
        maximum distance any tip can travel at each time step.

        * effect: Will label the tips who will annihilate (collide) with collision flag and also give it the address of the
        collision (the corresponding duct object that causes the collision as well as the position of the collision).

        * position of the collision: There are many choices but here we will use the closest point on the duct to the colliding tip

        * mechanism: We first create a KDtree of the whole structure. Then for each duct tip we find the points in the tree
        are in the r_c (collision radius) if the corresponding tip. Note all the points in the ball will cause annihilation of the
        tip. To clarify this, note that the annihilation between a tip cell and a duct cell is one of the following two cases:
        1) the duct cell is just created (in the last 2-3 steps) 2) the duct cell not just created.
        If the duct cell causing collision just created, then it must not be of the same branch or any of the sister branches. That
        is because there are always some duct that are just created by the moving tip (and thus are in the ball of collision) and
        also very shortly after branching the sister branches are very close to each other and we do not want this to be considered
        as a collision.
        For the second case, i.e. when the duct cell causing the collision is not recently created, then the collision happens
        without any special cases.
        """
        wholeNetwork = []
        pointIndexToDuct = []
        aliveDucts = []
        for duct in self.ductsList:
            wholeNetwork += duct.tailCurve
            pointIndexToDuct += [duct]*len(duct.tailCurve)
            if duct.tip.isAlive:
                aliveDucts.append(duct)
        numberOfAliveTips = len(aliveDucts)

        wholeNetwork = np.array(wholeNetwork)
        kdtree = KDTree(wholeNetwork)

        for duct in aliveDucts:
            nearbyIndices = kdtree.query_ball_point(duct.tip.position, self.r_c)
            for index in nearbyIndices:
                stopperDuct = pointIndexToDuct[index]
                stopperPoint = wholeNetwork[index]

                ## If the stopper duct is a sister duct, then at the first 10 steps we don't want to have
                # collision
                if stopperDuct.parent == duct.parent:
                    if stopperPoint.tolist() in stopperDuct.tailCurve[:10]:
                        continue

                ## if the stopper duct is the same duct or the parent duct then we don't want
                # to have a collision if the stopper Point is recently generated
                if stopperDuct == duct or stopperDuct == duct.parent:
                    if stopperPoint.tolist() in stopperDuct.tailCurve[-10:]:
                        continue


                # ## Check to see if the stopper is newly generated
                # if stopperPoint.tolist() in stopperDuct.tailCurve[-2:] : # the stopper is a newly generated point
                #     if stopperDuct.parent == duct.parent or stopperDuct == duct or stopperDuct == duct.parent:
                #         ## no annihilation happens
                #         continue

                duct.tip.collisionData = [stopperDuct, stopperPoint]
                break

        # plt.scatter(wholeNetwork[:,0],wholeNetwork[:,1])
        # plt.show()
