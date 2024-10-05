import numpy as np
import matplotlib.pyplot as plt
from agentBased.Network import Network



class Main:
    def __init__(self):
        self.network = Network()

    def run(self):
        for i in range(800):
            self.network.update()

    def plot(self):
        for duct in self.network.ductsList:
            coord = np.array(duct.tailCurve)
            # plt.scatter(coord[:,0],coord[:,1],s=2,c="k")
            plt.plot(coord[:,0],coord[:,1],color='steelblue')
            if duct.tip.isAlive:
                plt.plot(duct.tip.position[0],duct.tip.position[1],'o',color="orange",ms=3)
        plt.axis('equal')
        plt.show()




if __name__ == '__main__':
    main = Main()
    main.run()
    main.plot()
    # main.network.collisionDetection()



