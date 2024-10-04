import numpy as np
import matplotlib.pyplot as plt
from agentBased.Network import Network



class Main:
    def __init__(self):
        self.network = Network()

    def run(self):
        for i in range(100):
            self.network.update()

    def plot(self):
        for duct in self.network.ductsList:
            coord = np.array(duct.tailCurve)
            plt.plot(coord[:,0],coord[:,1])
        plt.axis('equal')
        plt.show()




if __name__ == '__main__':
    main = Main()
    main.run()
    main.plot()



