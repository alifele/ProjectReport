import numpy as np
import matplotlib.pyplot as plt
from agentBased.Network import Network
import networkx as nx
import sys

# Increase the recursion limit
sys.setrecursionlimit(10000)


class Main:
    def __init__(self):
        self.network = Network()
        self.G = None

    def run(self):
        for i in range(200):
            self.network.update()
        self.turnToGraph()
        self.G.remove_edges_from(nx.selfloop_edges(self.G))

    def plot(self):
        plt.figure(figsize=(12,12))
        for duct in self.network.ductsList:
            coord = np.array(duct.tailCurve)
            # plt.scatter(coord[:,0],coord[:,1],s=2,c="k")
            plt.plot(coord[:,0],coord[:,1],color='steelblue',alpha=0.2)
            if duct.tip.isAlive:
                pass
                # plt.plot(duct.tip.position[0],duct.tip.position[1],'o',color="orange",ms=3)

        # pos = {node: node for node in self.G.nodes()}
        # nx.draw(self.G, pos, node_size=5, node_color="lightblue", edge_color="red", font_size=8)

        cycles = nx.cycle_basis(self.G)
        cmap = plt.get_cmap('tab20')

        for idx, cycle in enumerate(cycles):
            color = cmap(idx % 20)
            for i in range(len(cycle) - 1):
                edge = main.G.get_edge_data(cycle[i], cycle[i + 1])
                curve = np.array(edge.get("curve"))
                plt.plot(curve[:, 0], curve[:, 1],color=color)

            i += 1
            edge = main.G.get_edge_data(cycle[i], cycle[0])
            curve = np.array(edge.get("curve"))
            plt.plot(curve[:, 0], curve[:, 1],color=color)

        plt.axis('equal')
        plt.show()

    def turnToGraph(self):
        self.G = nx.Graph()
        # Keep track of processed ducts to avoid duplicate recursion
        processed = set()
        for duct in self.network.ductsList:
            self.process_duct(duct, self.G, processed)

    def process_duct(self, duct, G, processed):
        """
        Process the duct and its connections recursively, ensuring no redundant processing of already visited ducts.
        """
        # Get the start point of the duct
        start_point = tuple(duct.tailCurve[0])

        if duct in processed:
            return
        processed.add(duct)

        # Process duct collisions (if any)
        if len(duct.collisionList) > 0:
            # Sort the collision list based on proximity to start_point
            collision_positions = np.array([tip.position for tip in duct.collisionList])
            pt_np = np.array(start_point)
            distances = np.sum((collision_positions - pt_np) ** 2, axis=1)
            sorted_indices = np.argsort(distances)
            sorted_collision_list = [duct.collisionList[i] for i in sorted_indices]

            for colliding_tip in sorted_collision_list:
                collision_point = tuple(colliding_tip.position)
                if collision_point not in G.nodes:
                    G.add_node(collision_point)

                # Add the edge between the start point and collision point
                curveStartIndex = duct.tailCurve.index(list(start_point))
                curveEndIndex = duct.tailCurve.index(list(collision_point))
                G.add_edge(start_point, collision_point,curve=duct.tailCurve[curveStartIndex:curveEndIndex+1])

                # Recursively process the colliding duct
                self.process_duct(colliding_tip.duct, G, processed)

                start_point = collision_point  # Update for next branch

        # Add the final end point of the duct
        end_point = tuple(duct.tailCurve[-1])
        if end_point not in G.nodes:
            G.add_node(end_point)

        # Create an edge between start and end if not already connected
        if not G.has_edge(start_point, end_point):
            curveStartIndex = duct.tailCurve.index(list(start_point))
            curveEndIndex = duct.tailCurve.index(list(end_point))
            G.add_edge(start_point, end_point,curve=duct.tailCurve[curveStartIndex:curveEndIndex+1])

    def turnToGraph_OLD(self):
        self.G = nx.Graph()
        for duct in self.network.ductsList:
            self.add_duct_to_graph(duct, self.G)

    def add_duct_to_graph_OLD(self, duct, G):
        # Add the start point of the duct as a node
        startPoint = tuple(duct.tailCurve[0])  # Assuming the first point in the curve is the start
        G.add_node(startPoint)

        # Explore collisions and add edges between start points of ducts
        if len(duct.collisionList) == 0:
            endPoint = tuple(duct.tailCurve[-1])
            G.add_edge(startPoint, endPoint)
        else:
            positions = np.array(
                [obj.position for obj in duct.collisionList])  # Shape (n, 2), where n is the number of objects
            pt_np = np.array(duct.tailCurve[0])  # Convert pt to a NumPy array
            distances = np.sum((positions - pt_np) ** 2, axis=1)
            sorted_indices = np.argsort(distances)
            sorted_collisionList = [duct.collisionList[i] for i in sorted_indices]

            # for colliding_tip in duct.collisionList:
            for colliding_tip in sorted_collisionList:
                # Get the start point of the colliding duct
                collisionPoint = tuple(colliding_tip.position)
                G.add_node(collisionPoint)
                G.add_edge(startPoint, collisionPoint)
                startPoint = collisionPoint
                self.add_duct_to_graph(colliding_tip.duct, G)

            endPoint = tuple(duct.tailCurve[-1])
            G.add_node(endPoint)
            G.add_edge(startPoint, endPoint)


if __name__ == '__main__':
    main = Main()
    main.run()
    main.plot()
    # main.network.collisionDetection()



