import numpy as np

from agentBased.utils import calculate_curve_length
import networkx as nx
from scipy.sparse import linalg as spla



class Remodel(object):
    def __init__(self, G):
        self.sinkNodes = None
        self.sourceNodes = None
        self.N_source = 1
        self.N_sink = 200
        self.G = G
        self.initializeGraph(self.G)
        self.V = self.calculateVolume(self.G)

        self.dt = 0.1   ## This dt is for the remodeling
        self.beta = 1
        self.mu = 1




    def remodel(self):
        for i in range(60):
            self.G = self.calculateFlow(self.G)
            self.G = self.updateConductance(self.G)


    def initializeGraph(self, G):
        # Step 1: Iterate over the edges and add attributes C, Q, and length
        for u, v, data in G.edges(data=True):
            data['C'] = 1.0  # Add attribute C with value 1.0
            data['Q'] = 0.0  # Add attribute Q with value 0.0
            # Assume 'curve' is a list of points (e.g., np.array with shape (N, 2) or a list of coordinates)
            if 'curve' in data:
                data['length'] = calculate_curve_length(data['curve'])  # Set 'length' to be the length of the curve
            else:
                data['length'] = 0  # If no 'curve', set length to 0 or any other default value

        # Step 2: Iterate over the nodes and assign attributes p=0.0 and q=0.0
        for node in G.nodes():
            G.nodes[node]['p'] = 0.0  # Add or set p attribute
            # G.nodes[node]['q'] = 0.0  # Add or set q attribute
            G.nodes[node]["s"] = 0.0

        # self.randomSourceSink(G)
        self.annulusSourceSink(G)




    def calculateVolume(self,G):
        total_volume = 0
        for u, v, data in G.edges(data=True):
            L = data["length"]
            C = data['C']
            total_volume += L * np.sqrt(C)
        return total_volume


    def g(self, Q):
        return np.abs(Q) ** (2 / 3)


    def calculateFlow(self, G):
        L = nx.laplacian_matrix(G, weight="C")  # Convert to a dense NumPy array if needed

        # Perturb the diagonal entry of the Laplacian matrix for one of the sink nodes
        sink_node = list(G.nodes).index(self.sinkNodes[np.random.randint(self.N_sink)]) # Choose one sink node for perturbation
        L = L.tolil()  # Convert to LIL format to modify the sparse matrix
        L[sink_node, sink_node] = L[sink_node, sink_node] + 1e-5  # Add a small constant to the diagonal
        L = L.tocsr()  # Convert back to CSR format for efficient solving

        S = np.array([G.nodes[n]['s'] for n in G.nodes()])
        P = spla.spsolve(L, S)

        # Assign pressures back to the graph nodes
        for i, n in enumerate(G.nodes()):
            G.nodes[n]['p'] = P[i]

        for u, v, data in G.edges(data=True):
            C = G.edges[u, v]["C"]
            # pos_u = np.array(G.nodes[u]['pos'])
            # pos_v = np.array(G.nodes[v]['pos'])
            # L = np.linalg.norm(pos_u - pos_v)
            L = data["length"]
            G.edges[u, v]["Q"] = C * (G.nodes[u]["p"] - G.nodes[v]["p"]) / L

        return G


    def updateConductance(self, G):
        C = np.array([G.edges[edge]['C'] for edge in G.edges()])
        Q = np.array([G.edges[edge]['Q'] for edge in G.edges()])
        A = 0
        for u, v, data in G.edges(data=True):
            # pos_u = np.array(G.nodes[u]['pos'])
            # pos_v = np.array(G.nodes[v]['pos'])
            # L = np.linalg.norm(pos_u - pos_v)
            L = data["length"]
            A += L * self.g(G.edges[u, v]["Q"])

        C = (np.sqrt(C) + self.dt * (self.V / self.beta * (self.g(Q) / A) - self.mu * np.sqrt(C))) ** 2

        for edge, new_conductance in zip(G.edges(), C):
            G.edges[edge]['C'] = new_conductance

        return G


    def randomSourceSink(self, G):
        self.sinkNodes = [list(G.nodes())[i] for i in np.random.choice(len(G.nodes()), self.N_sink, replace=False)]
        # self.sourceNodes = [list(G.nodes())[i] for i in np.random.choice(len(G.nodes()), self.N_source, replace=False)]
        self.sourceNodes = [(0.0, 0.0)]

        for n in self.sourceNodes:
            G.nodes[n]["s"] = 1 / self.N_source
        for n in self.sinkNodes:
            G.nodes[n]["s"] = -1 / self.N_sink

        return G


    def annulusSourceSink(self, G):
        self.sourceNodes = [(0.0, 0.0)]

        nodeDistanceFromOriginList = []
        for node in G.nodes():
            nodeDistanceFromOriginList.append(np.linalg.norm(node))
        R_max = np.max(nodeDistanceFromOriginList)
        R_min = R_max / 2 ##

        candidateNodes = []
        for i, node in enumerate(G.nodes()):
            if R_min < nodeDistanceFromOriginList[i] < R_max :
                candidateNodes.append(node)

        indices = np.random.choice(len(candidateNodes), size=self.N_sink, replace=False)  # Choose n unique random indices
        self.sinkNodes =  [candidateNodes[i] for i in indices]


        for n in self.sourceNodes:
            G.nodes[n]["s"] = 1 / self.N_source
        for n in self.sinkNodes:
            G.nodes[n]["s"] = -1 / self.N_sink






