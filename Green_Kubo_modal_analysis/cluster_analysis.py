
# 2019-05-01  YJ
# given a neighbor list, e.g. [[1,2],[2,0], [0, 1]]  where 0 is connected to 1 and 2, 1 is connected to 0 and 2, 2 is connected to 0 and 1
# calculate the number of nth neighbor that is far from a central atom but can be reach by n step via a shortest path
# implement BSF search
# output a list of number of nth neighbor [] index means nth neighbor, 0th is the central atom itself.

# Python3 Program to print BFS traversal
# from a given source vertex. BFS(int s)
# traverses vertices reachable from s.


# This class represents a directed graph
# using adjacency list representation
# This code is adapted from Neelam Yadav

from collections import defaultdict
from basic_function import bondlen

class Graph:

    # Constructor
    def __init__(self):
        # default dictionary to store graph
        self.graph = defaultdict(list)

    # function to add an edge to graph
    def addEdge(self, u, v):
        self.graph[u].append(v)

    # add all edge from neighbor list
    def add_all_Edge(self, neighborlist):
        for i in neighborlist:
            for j in neighborlist[i]:
                self.graph[i].append(j)

    # Function to print a BFS of graph
    def BFS(self, s):
        # Mark all the vertices as not visited
        visited = [False] * (len(self.graph))

        # Create a queue for BFS
        queue = []

        # Mark the source node as
        # visited and enqueue it
        queue.append(s)
        visited[s] = True

        while queue:
            # Dequeue a vertex from
            # queue and print it
            s = queue.pop(0)
            # Get all adjacent vertices of the
            # dequeued vertex s. If a adjacent
            # has not been visited, then mark it
            # visited and enqueue it
            for i in self.graph[s]:
                if visited[i] == False:
                    queue.append(i)
                    visited[i] = True

    def gen_nth_neighbor(self, mysystem, frame, maximum_neighbor_order = 8):

        lx, ly, lz = mysystem[frame].L[0], mysystem[frame].L[1], mysystem[frame].L[2]
        res = [1.0]
        distance = [0.0]
        count_bond_to_next_order = [0.0]
        if bool(self.graph) is False:
            print("neighbor list is empty!")
            return res, distance
        for _ in range(1, maximum_neighbor_order + 1):
            res.append(0.0)
            distance.append(0.0)
            count_bond_to_next_order.append(0.0)

        visited = set()
        si = []
        # the following 4 lines just help debug
        for i in self.graph:
            si.append(i)
        si.sort()
        for n in si:
            # print("n=",n)
            visited.clear()
            # Mark all the vertices as not visited
            # Create a queue for BFS
            queue = []
            # Mark the source node as
            # visited and enqueue it
            queue.append(n)
            visited.add(n)
            order = 1
            size = 1

            next_level_node = set()

            while queue and order <= maximum_neighbor_order:
                # Dequeue a vertex from
                # queue and print it
                next_level_node.clear()
                count_bond = 0
                count = 0
                dist = 0.0
                for _ in range(size):
                    s = queue.pop(0)
                    for i in self.graph[s]:
                        if i in next_level_node:
                            count_bond += 1
                        if not i in visited:
                            count_bond += 1
                            queue.append(i)
                            visited.add(i)
                            next_level_node.add(i)
                            count += 1
                            dist += bondlen(mysystem[frame].myatom[n], mysystem[frame].myatom[i], lx, ly, lz)
                res[order] += count
                count_bond_to_next_order[order] += count_bond
                if count > 0:
                    distance[order] += dist / count
                size, order = count, order + 1
        print("there are:", len(self.graph), "central atoms")
        for n in range(1, maximum_neighbor_order + 1):
            res[n] /= len(self.graph)
            distance[n] /= len(self.graph)
            count_bond_to_next_order[n] /= len(self.graph)

        return res, distance, count_bond_to_next_order


    # 2020-04-23 get coordination sequence (CSQ) according to paper by Meier and Moeck, 1979, J. solid state chemistry.
    # would return all possible CSQ, their N2, N5, N6 values and their frequency which can be visualized using color
    # CSQ store the coordination sequence from each Si atom in format: tuple (N1,N2,N3,N4,N5,N6) as key
    # and [N1, N2, N5, N6, Total(nodes within the cluster, N6 shared by others so divided by two), frequency]
    def gen_nth_neighbor_CSQ(self, mysystem, CSQ, frame, record_atom_density, maximum_neighbor_order = 6):

        assert maximum_neighbor_order == 6 # make sure ==6, otherwise need to modify code, because the 6th element will be read in the CSQ element list

        lx, ly, lz = mysystem[frame].L[0], mysystem[frame].L[1], mysystem[frame].L[2]
        if bool(self.graph) is False:
            print("neighbor list is empty!")
            return

        visited = set()
        si = []
        # the following 4 lines just help debug
        for i in self.graph:
            si.append(i)
        si.sort()

        #CSQ = defaultdict(list)
        single_CSQ = [] # for each atom
        for n in si:
            # print("n=",n)
            single_CSQ = [0] * maximum_neighbor_order
            visited.clear()
            # Mark all the vertices as not visited
            # Create a queue for BFS
            queue = []
            # Mark the source node as
            # visited and enqueue it
            queue.append(n)
            visited.add(n)
            order = 1
            size = 1

            next_level_node = set()
            ave_N6_distance = 0
            while queue and order <= maximum_neighbor_order:
                # Dequeue a vertex from
                # queue and print it
                next_level_node.clear()
                count = 0
                dist = 0.0
                for _ in range(size):
                    s = queue.pop(0)
                    for i in self.graph[s]:
                        if not i in visited:
                            queue.append(i)
                            visited.add(i)
                            next_level_node.add(i)
                            count += 1
                            dist += bondlen(mysystem[frame].myatom[n], mysystem[frame].myatom[i], lx, ly, lz)
                single_CSQ[order - 1] += count
                if order == maximum_neighbor_order:
                    ave_N6_distance = dist / count
                size, order = count, order + 1
            assert single_CSQ[5] > 0
            vol_N6 = 4.0/3.0 * 3.1415926 * ave_N6_distance ** 3
            record_atom_density[n] = [single_CSQ[0], (1.0 + sum(single_CSQ) - 0.50 * single_CSQ[5]) / vol_N6,
                                      ave_N6_distance, single_CSQ[2] + single_CSQ[3], single_CSQ[0] + single_CSQ[1]]
            CSQ_code = tuple(single_CSQ)
            #CSQ_code2 = tuple(single_CSQ[0:2])
            #CSQ_code4 = tuple(single_CSQ[0:2])
            # treat cluster as the same no matter their difference in radius which will be got from average of all distance
            if CSQ_code in CSQ:
                #  add one more count/frequency to the same CSQ species
                CSQ[CSQ_code][-2] += 1
                CSQ[CSQ_code][-1] += ave_N6_distance
            else:
                # store N1, N2, N5, N6, (N1 + N2 +... + N6)
                total_N1_N2 = single_CSQ[0] + single_CSQ[1]
                total_N3_N4 = single_CSQ[2] + single_CSQ[3]
                total_N1_N6 = 1.0 + sum(single_CSQ) - 0.50 * single_CSQ[5]
                CSQ[CSQ_code] = [single_CSQ[0], single_CSQ[1], single_CSQ[4], single_CSQ[5], total_N1_N6, total_N1_N2, total_N3_N4, 1, ave_N6_distance]
        print("there are:", len(self.graph), "central atoms")
        #for n in CSQ.items():
        #    CSQ[n][4] /= len(self.graph)
        return
