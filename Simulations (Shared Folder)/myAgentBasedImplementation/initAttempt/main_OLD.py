import matplotlib.pyplot as plt
from numpy import linalg as LA
from pylab import *
import networkx as nx
from scipy.spatial import KDTree




class BARW:
    def __init__(self,tmax):
        self.G = nx.Graph()
        self.coordinates = None
        self.evolve = None
        self.diff = None
        self.angle_list = None
        self.angle = None
        self.angle_init = None
        self.node = None
        self.current_edge = {}
        self.active_tips = {}  # Track current edges for each active tip
        self.branch_counter = 0  # To uniquely identify branches
        self.coordsTemp = None

        self.branchID_TO_StartNode = {} # Give the branchID as the key, it will give you the StartNode
        self.branchID_TO_endNode = {}


        ## Parameters
        self.Lx = 200
        self.Lz = 300
        self.tmax = tmax

        self.lstep = 1                  # Elongation length (step size)
        self.min_branch = pi / 10       # minimum branching angle
        self.cutoff = 1.5 * self.lstep  # cutoff length for annihilation via tip-node contact
        self.radavoid = 10 * self.lstep # Radius and strength of attraction/repulsion (if fav<0, then attraction)



    def initialize(self):
        self.node = np.array([[0, self.Lz / 2, 0, 1]])      # [x,y,genID, branchID] for the tip cells only
        self.angle_init = [0.0]
        self.angle = np.array(self.angle_init)              # angle for the tip cells only
        self.angle_list = np.array([[0, 1]])                # [angle, branchID] TODO: the second entry should be 1?
        self.diff = 0
        self.evolve = np.array([len(self.node)])
        self.coordinates = np.array([self.node[0]])

        # initial_point = (self.node[0][0], self.node[0][1])
        initial_point = tuple(self.node[0][:2])
        self.G.add_node(initial_point)  # Add the initial point as the first node
        self.branchID_TO_StartNode[1] = initial_point
        # self.current_edge[0] = initial_point  # Track the current starting point of the first tip's branch

    def distance(self, vect, tip):
        diff = np.add(vect[:2], -tip[:2])
        distance = LA.norm(diff)
        return distance


    def prob_list(self, x, f, pb):
        # bias for forward steps
        alpha = np.clip((1 - f * np.sin(x)) / 2, 0, 1)
        # bias for backward steps
        beta = np.clip((1 + f * np.sin(x)) / 2, 0, 1)
        # list of stepping probabilities:
        problist = np.array([(1 - pb) * alpha, (1 - pb) * beta, pb * alpha, pb * beta])
        # 'cumulative' probabilities of the list
        cum_prob = [0] + [sum(problist[:j + 1]) for j in range(len(problist))]
        probs = {'problist': problist, 'cumulative': np.array(cum_prob)}
        return probs

    def self_avoidance(self, fav):
        """
        Implements the self-avoidance rule using KDTree for efficiency.
        Avoidance potential is applied only if fav != 0.
        """
        if fav == 0:
            return

        passive_coords = self.coordinates[:, :2]  # Take only (x, y) coordinates
        if len(passive_coords) > 0:
            kdtree = KDTree(passive_coords)

        for j, tip in enumerate(self.node):
            # Query the KDTree for points within the avoidance radius from the tip
            nearby_indices = kdtree.query_ball_point(tip[:2], self.radavoid)
            dist_sum = np.zeros(2)

            for idx in nearby_indices:
                # Ignore distances to parent, sister branches, and nodes within the same duct
                if tip[-2] == self.coordinates[idx, -1] or tip[-1] == self.coordinates[idx, -1] or tip[-2] == \
                        self.coordinates[idx, -2]:
                    continue

                distance = self.distance(self.coordinates[idx], tip)
                if 0 < distance <= self.radavoid:
                    # Sum of displacement vectors within the avoidance potential radius
                    dist_sum += self.coordinates[idx, :2] - tip[:2]

            # Normalize and apply the avoidance potential
            norm_dis = LA.norm(dist_sum)
            if norm_dis > 0:
                displace = fav * (dist_sum / norm_dis)
                self.node[j][:2] += displace  # Update node position with displacement

    def annihilation(self):
        """
        Implements the annihilation rule using KDTree for efficiency.
        Active tips are removed if they are within the cutoff distance of passive nodes.
        """
        newNodes = []

        passive_coords = self.coordinates[:, :2]  # Take only (x, y) coordinates
        if len(passive_coords) > 0:
            kdtree = KDTree(passive_coords)

        skipp = 0
        # List of recent nodes generated in the last 2 time steps
        checklist = [list(item) for item in self.coordinates[-int(np.sum(self.evolve[-2:])):]]

        for j in range(len(self.node)):  # Loop over the tip cells
            tip = np.array(self.node[j - skipp])
            if len(self.node) == 0:
                break

            # Query the KDTree for points within the cutoff distance from the tip
            nearby_indices = kdtree.query_ball_point(tip[:2], self.cutoff)

            once = False
            for idx in nearby_indices:
                nearby_node = passive_coords[idx]
                radius = self.distance(nearby_node, tip[:2])

                # Check if the nearby node belongs to the same duct, parent, or sister branch
                if list(self.coordinates[idx]) in checklist:
                    if self.coordinates[idx, -1] != tip[-1] and self.coordinates[idx, -1] != tip[-2]: ## first: not the same branch, second: not parent
                        if 0 < radius < self.cutoff:
                            self.node = np.delete(self.node, j - skipp, 0)
                            self.angle = np.delete(self.angle, j - skipp, 0)
                            once = True
                            skipp += 1
                elif 0 < radius < self.cutoff:
                    self.node = np.delete(self.node, j - skipp, 0)
                    self.angle = np.delete(self.angle, j - skipp, 0)
                    once = True
                    skipp += 1

                if once:
                    # self.coordinates[idx][:2] = nearby_node ## Make the annator point the same as the annated tip node. This
                    # makes sure that the the curves are joined

                    # self.node[j - skipp][:2] = nearby_node[:2] ## Make sure that the tip is considered the same as the
                    # annator nearby point. This will make the curves to join
                    coord_temp = np.append(self.coordinates, np.array(self.node), axis=0)
                    topnr = max(coord_temp[:, -1]) ## To get the highest branchID

                    annatedEndNode = nearby_node[:2]
                    annatedBranchID = tip[3]
                    annatedGenID = tip[2]
                    annatedStartNode = self.branchID_TO_StartNode[annatedBranchID]
                    self.branchID_TO_endNode[annatedBranchID] = annatedEndNode
                    self.G.add_node(tuple(annatedEndNode))
                    self.G.add_edge(tuple(annatedStartNode), tuple(annatedEndNode))
                    annated_coordEntry = [annatedEndNode[0], annatedEndNode[1], annatedGenID, annatedBranchID]
                    # self.coordinates = np.append(self.coordinates, [annated_coordEntry], axis=0)
                    self.coordsTemp = np.append(self.coordsTemp, [annated_coordEntry], axis=0)

                    ## TODO: Need to add angle entries as well

                    annatorBranchID = self.coordinates[idx][3]
                    annatorGenID = self.coordinates[idx][2]
                    annatorStartNode = self.branchID_TO_StartNode[annatorBranchID]
                    if annatorBranchID in self.branchID_TO_endNode.keys(): ## The tip cell of annator is dead
                        annatorEndNode = self.branchID_TO_endNode[annatorBranchID]
                        self.G.remove_edge(tuple(annatorStartNode), tuple(annatorEndNode))
                        # self.G.remove_nodes_from([tuple(annatorStartNode), tuple(annatorEndNode)])
                        del self.branchID_TO_endNode[annatorBranchID]
                        del self.branchID_TO_StartNode[annatorBranchID]

                    sec1StartNode = annatorStartNode
                    sec1EndNode = annatedEndNode
                    sec1BranchID = annatorBranchID
                    sec1GenID = annatorGenID
                    sec1_coordEntry = [sec1EndNode[0],sec1EndNode[1],sec1GenID, sec1BranchID]
                    # self.coordinates = np.insert(self.coordinates, idx, sec1_coordEntry, axis=0)
                    self.coordsTemp = np.insert(self.coordsTemp, idx, sec1_coordEntry, axis=0)
                    self.branchID_TO_StartNode[sec1BranchID] = sec1StartNode
                    self.branchID_TO_endNode[sec1BranchID] = sec1EndNode
                    self.G.add_nodes_from([tuple(sec1StartNode), tuple(sec1EndNode)])
                    self.G.add_edge(tuple(sec1StartNode), tuple(sec1EndNode))
                    # self.coordinates[:idx+1, 3] = np.where(self.coordinates[:idx+1, 3] == annatorBranchID, sec1BranchID, self.coordinates[:idx+1, 3])
                    self.coordsTemp[:idx+1, 3] = np.where(self.coordsTemp[:idx+1, 3] == annatorBranchID, sec1BranchID, self.coordsTemp[:idx+1, 3])


                    if np.sum(self.coordsTemp[idx:,3]==annatorBranchID) > 0 : ## Then there is a second section as well
                        sec2StartNode = annatedEndNode
                        sec2BranchID = topnr + 1
                        # self.coordinates[idx+1:, 3] = np.where(self.coordinates[idx+1:, 3] == annatorBranchID,
                        #                                      sec2BranchID, self.coordinates[idx+1:, 3])
                        self.coordsTemp[idx+1:, 3] = np.where(self.coordsTemp[idx+1:, 3] == annatorBranchID,
                                                             sec2BranchID, self.coordsTemp[idx+1:, 3])
                        self.branchID_TO_StartNode[sec2BranchID] = tuple(sec2StartNode)

                        if annatorBranchID in self.node[:,3]: ## This means that the tip of annator is still alive
                            # self.node[self.node[:,3] == annatorBranchID][0][-1] = sec2BranchID
                            myNode = self.node[self.node[:,3]==annatorBranchID][0]
                            # self.node[self.node[:,3] == annatorBranchID, -1] = sec2BranchID
                            newNodes.append([myNode[0], myNode[1], myNode[2], sec2BranchID])
                        else:
                            sec2EndNode = self.coordsTemp[idx:][self.coordsTemp[idx:,3] == annatorBranchID][-1][:2]
                            self.G.add_nodes_from([tuple(sec2StartNode), tuple(sec2EndNode)])
                            self.G.add_edge(tuple(sec2StartNode), tuple(sec2EndNode))
                            # self.branchID_TO_StartNode[sec2BranchID] = tuple(sec2StartNode)
                            self.branchID_TO_endNode[sec2BranchID] = tuple(sec2EndNode)

                    break
        for node in newNodes:
            # self.node[self.node[:,:2] == node[:2], -1] = node[3]
            row_index = np.where((self.node[:, 0] == node[0]) & (self.node[:, 1] == node[1]))[0]
            # If a matching row is found, set the last element to node[3]
            if row_index.size > 0:
                self.node[row_index[0], -1] = node[3]





    def tissue1(self, prob_branch, fav, fchem):
        ## By the following two lists we want to make sure that at the branching point or at the splitting point of annator
        # the branching point or the annator point acts as the common end/start point for the parent and daughter branches
        coordEntryList = []
        angleEntryList = []
        branchedFlag = False
        self.coordsTemp = np.copy(self.coordinates)
        skip = 0
        rnr = rand(len(self.angle))
        cumprob = np.array(
            [self.prob_list(np.radians(item), fchem, prob_branch)['cumulative'] for item in self.angle_list[-(len(self.angle)):, 0]])
        index_rnr = np.array([np.where(cumprob[j] > rnr[j])[0][0] - 1 for j in range(len(rnr))])

        node_temp = self.coordinates[-int(self.evolve[-1]):]

        for j in range(len(index_rnr)):  # Loop over all active tip cells.
            coord_temp = np.append(self.coordinates, np.array(self.node), axis=0)
            topnr = max(coord_temp[:, -1])

            if index_rnr[j] == 0:
                ang_elong = np.random.uniform(0, self.min_branch)
                self.angle[j + skip] += ang_elong
                self.node[j + skip] = [self.node[j + skip][0] + self.lstep * cos(self.angle[j + skip]), \
                                  self.node[j + skip][1] + self.lstep * sin(self.angle[j + skip]), self.node[j + skip][2], self.node[j + skip][3]]
            elif index_rnr[j] == 1:
                ang_elong = np.random.uniform(0, self.min_branch)
                self.angle[j + skip] -= ang_elong
                self.node[j + skip] = [self.node[j + skip][0] + self.lstep * cos(self.angle[j + skip]), \
                                  self.node[j + skip][1] + self.lstep * sin(self.angle[j + skip]), self.node[j + skip][2], self.node[j + skip][3]]
            elif index_rnr[j] >= 2:
                branchedFlag = True
                ang_branch1 = np.random.uniform(self.min_branch, pi / 2)
                ang_branch2 = np.random.uniform(self.min_branch, pi / 2)

                ## Make sure to add the EndNode of the Mother and draw the edge
                motherBranchID = self.node[j + skip][3]
                motherEndNode = tuple(self.node[j + skip][:2])
                motherStartNode = self.branchID_TO_StartNode[motherBranchID]
                motherAngle = self.angle[j+skip]
                self.branchID_TO_endNode[motherBranchID] = motherEndNode
                self.G.add_node(motherEndNode)
                self.G.add_edge(motherStartNode, motherEndNode)

                self.angle = np.insert(self.angle, j + skip + 1, self.angle[j + skip] + ang_branch1)
                self.node = np.insert(self.node, j + skip + 1, [self.node[j + skip][0] + self.lstep * cos(self.angle[j + skip + 1]), \
                                                      self.node[j + skip][1] + self.lstep * sin(self.angle[j + skip + 1]),
                                                      self.node[j + skip][3], topnr + 2], axis=0)
                self.angle[j + skip] = self.angle[j + skip] - ang_branch2
                self.node[j + skip] = [self.node[j + skip][0] + self.lstep * cos(self.angle[j + skip]), \
                                  self.node[j + skip][1] + self.lstep * sin(self.angle[j + skip]), self.node[j + skip][3], topnr + 1]


                daughter1_StartNode = motherEndNode
                daughter2_StartNode = motherEndNode
                daughter1_BranchID = self.node[j + skip][3]
                daughter2_BranchID = self.node[j + skip + 1][3]
                daughter1_GenID = self.node[j + skip][2]
                daughter2_GenID = self.node[j + skip + 1][2]
                daughter1_angle = [np.degrees(motherAngle),daughter1_BranchID]
                daughter2_angle = [np.degrees(motherAngle),daughter2_BranchID]
                self.branchID_TO_StartNode[daughter1_BranchID] = daughter1_StartNode
                self.branchID_TO_StartNode[daughter2_BranchID] = daughter2_StartNode
                self.G.add_node(daughter1_StartNode)
                self.G.add_node(daughter2_StartNode)

                ## We want to make sure we submit coordinate entries of daughter 1 and daughter 2 with their start node
                # position as the mother end node (not the updated corresponding tips above)
                daughter1_coordEntry = [daughter1_StartNode[0], daughter1_StartNode[1],daughter1_GenID, daughter1_BranchID]
                daughter2_coordEntry = [daughter2_StartNode[0], daughter2_StartNode[1],daughter2_GenID, daughter2_BranchID]
                self.coordsTemp = np.append(self.coordsTemp, np.array([daughter1_coordEntry, daughter2_coordEntry]), axis=0)
                # coordEntryList.append(daughter1_coordEntry)
                # coordEntryList.append(daughter2_coordEntry)
                # angleEntryList.append(daughter1_angle)
                # angleEntryList.append(daughter2_angle)

                # self.coordinates = np.append(self.coordinates, [daughter1_coordEntry],axis=0)
                # self.coordinates = np.append(self.coordinates, [daughter2_coordEntry],axis=0)


                # self.record_branching(self.node[j + skip], self.node[j + skip + 1], topnr)  # Record branching
                skip += 1

        self.angle = (self.angle + pi) % (2 * pi) - pi
        self.self_avoidance(fav)
        self.annihilation()
        self.evolve = np.append(self.evolve, len(self.node))
        self.angle = (self.angle + pi) % (2 * pi) - pi
        # if branchedFlag == True:
        #     self.coordinates = np.append(self.coordinates, np.array(coordEntryList),axis=0)
        #     self.angle_list = np.append(self.angle_list, np.array(angleEntryList),axis=0)

        self.angle_list = np.append(self.angle_list, np.column_stack((np.degrees(self.angle), self.node[:, -1])), axis=0)
        # self.coordsTemp = np.append(self.coordsTemp, np.array(self.node), axis=0)
        # self.coordinates = np.copy(self.coordsTemp)
        self.coordinates = np.append(self.coordinates, np.array(self.node), axis=0)



    def record_branching(self, parent_tip, daughter_tip, branch_id):
        """ Records a branching event by creating new nodes and edges in the graph """
        parent_tuple = tuple(parent_tip[:2])
        daughter_tuple = tuple(daughter_tip[:2])
        if parent_tuple not in self.G:
            self.G.add_node(parent_tuple)
        self.G.add_node(daughter_tuple)
        self.G.add_edge(parent_tuple, daughter_tuple)

    def remove_arcs(self):
        # Identify arcs (self-loops)
        arcs = list(nx.selfloop_edges(self.G))
        # Remove the arcs from the graph
        self.G.remove_edges_from(arcs)

    def tissue1_OLD(self, prob_branch, fav, fchem):
        '''
        node: contains the [x,y,genID,branchID] of the tips only
        angle: contains the [\phi, branchID] of the tips only
        coordinates: contains the [x,y,genID,branchID] of every point of the branches
        angle_list: contains the [\phi, branchID] of every points of the branches
        '''
        skip = 0
        # draw random numbers to decide on next jumps:
        rnr = rand(len(self.angle))
        # determine the cumulative distribution of the stepping probabilities:
        # (use angle values from the local angle list to impose parallel field as guidance)
        cumprob = np.array(
            [self.prob_list(np.radians(item), fchem, prob_branch)['cumulative'] for item in self.angle_list[-(len(self.angle)):, 0]])
        # determine the first entry of cumprob that is larger than rnr, and take the entry before that:
        index_rnr = np.array([np.where(cumprob[j] > rnr[j])[0][0] - 1 for j in range(len(rnr))])

        node_temp = self.coordinates[-int(self.evolve[-1]):]

        for j in range(len(index_rnr)):  # This loops over all the active tip cells.
            coord_temp = np.append(self.coordinates, np.array(self.node), axis=0)
            # highest tree label of existing nodes
            topnr = max(coord_temp[:, -1])

            # elongation happens for the first two entries of the cumulative dist:
            if index_rnr[j] == 0:
                # determine a random elongation angle
                # change the angle and coordinates of active tips
                ang_elong = np.random.uniform(0, self.min_branch)
                # the first entry of cumulative distribution cooresponds to a forward step, i.e. increase the local angle!
                self.angle[j + skip] += ang_elong
                self.node[j + skip] = [self.node[j + skip][0] + self.lstep * cos(self.angle[j + skip]), \
                                  self.node[j + skip][1] + self.lstep * sin(self.angle[j + skip]), self.node[j + skip][2], self.node[j + skip][3]]
            elif index_rnr[j] == 1:
                # determine a random elongation angle
                # change the angle and coordinates of active tips
                ang_elong = np.random.uniform(0, self.min_branch)
                # the second entry of cumulative distribution cooresponds to a backward step, i.e. decrease the local angle!
                self.angle[j + skip] -= ang_elong
                self.node[j + skip] = [self.node[j + skip][0] + self.lstep * cos(self.angle[j + skip]), \
                                  self.node[j + skip][1] + self.lstep * sin(self.angle[j + skip]), self.node[j + skip][2], self.node[j + skip][3]]

                # branching happens for the last two entries of the cumulative dist:
            ## The following code is when we did not record any graph data strcure.
            elif index_rnr[j] >= 2:
                # determine two random angles between pi/10 and pi/2:
                ang_branch1 = np.random.uniform(self.min_branch, pi / 2)
                ang_branch2 = np.random.uniform(self.min_branch, pi / 2)
                # add a new branch changing the coordinates with the random angle ang_branch1:
                self.angle = np.insert(self.angle, j + skip + 1, self.angle[j + skip] + ang_branch1)
                self.node = np.insert(self.node, j + skip + 1, [self.node[j + skip][0] + self.lstep * cos(self.angle[j + skip + 1]), \
                                                      self.node[j + skip][1] + self.lstep * sin(self.angle[j + skip + 1]),
                                                      self.node[j + skip][3], topnr + 2], axis=0)
                # change the angle and coordinates of the remaining branch with the random angle ang_branch2:
                self.angle[j + skip] = self.angle[j + skip] - ang_branch2
                self.node[j + skip] = [self.node[j + skip][0] + self.lstep * cos(self.angle[j + skip]), \
                                  self.node[j + skip][1] + self.lstep * sin(self.angle[j + skip]), self.node[j + skip][3], topnr + 1]
                skip += 1

        # Call the self-avoidance function after node updates
        self.angle = (self.angle + pi) % (2 * pi) - pi
        self.self_avoidance(fav)
        self.annihilation()

        # save the length of node vector (to track node evolution over time)
        self.evolve = np.append(self.evolve, len(self.node))

        # set angle values to be within [-pi,pi]
        self.angle = (self.angle + pi) % (2 * pi) - pi
        # save the angles of nodes [in degrees!] including the generation number
        self.angle_list = np.append(self.angle_list, np.column_stack((np.degrees(self.angle), self.node[:, -1])), axis=0)
        # save the coordinates of all nodes
        self.coordinates = np.append(self.coordinates, np.array(self.node), axis=0)


    def plotter_basic(self):
        # Plot simulated network

        fig, ax = plt.subplots(figsize=(4, 4))
        ms = 1.5

        # this function can be used to plot the network until a certain time point

        def step(till, evolve):
            step = np.sum(evolve[:till])
            return int(step)

        # choose time1 to plot until a certain timepoint (time1=tmax for the complete network)
        time1 = self.tmax - 1
        time2 = time1 + 1

        x = [item[0] for item in self.coordinates[:step(time2, self.evolve)]]
        y = [item[1] for item in self.coordinates[:step(time2, self.evolve)]]

        ax.plot(x, y, 'o', color='steelblue', markersize=ms)
        ax.plot(x[0], y[0], 'x', color='firebrick', markersize=8)

        # plot active tips with different color
        ax.plot(x[step(time1, self.evolve):step(time2, self.evolve)], y[step(time1, self.evolve):step(time2, self.evolve)], 'o', color='C1',
                markersize=ms + 0.5)

        plt.savefig("pythonGen.png",dpi=300)

    def plot_graph(self):
        """
        Plot the graph structure generated by the BARW model.
        Each node represents either the initial point, a branching point, or an annihilation point.
        Each edge represents a duct (inactive tips) between nodes.
        """
        pos = {node: (node[0], node[1]) for node in self.G.nodes()}  # Position nodes based on their coordinates
        plt.figure(figsize=(8, 8))

        # Plot nodes
        nx.draw_networkx_nodes(self.G, pos, node_size=50, node_color='blue')

        # Plot edges
        nx.draw_networkx_edges(self.G, pos, edgelist=self.G.edges(), edge_color='black')

        # Set plot details
        plt.title("Graph Structure of the BARW Model")
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()
