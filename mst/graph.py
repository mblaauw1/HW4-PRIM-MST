import numpy as np
from typing import Union
import pytest

class Graph:

    def __init__(self, adjacency_mat: Union[np.ndarray, str]):
        """
    
        Unlike the BFS assignment, this Graph class takes an adjacency matrix as input. `adjacency_mat` 
        can either be a 2D numpy array of floats or a path to a CSV file containing a 2D numpy array of floats.

        In this project, we will assume `adjacency_mat` corresponds to the adjacency matrix of an undirected graph.
    
        """
        if type(adjacency_mat) == str:
            self.adj_mat = self._load_adjacency_matrix_from_csv(adjacency_mat)
        elif type(adjacency_mat) == np.ndarray:
            self.adj_mat = adjacency_mat
   
        
        else: 
            raise TypeError('Input must be a valid path or an adjacency matrix')
        self.mst = []

    def _load_adjacency_matrix_from_csv(self, path: str) -> np.ndarray:
        with open(path) as f:
            return np.loadtxt(f, delimiter=',')
    
   # def _load_adjacency_matrix_from_csv(self, path: str) -> np.ndarray:
    #    with open(path) as f:
     #       lines = f.readlines()

    # Assuming the matrix is rectangular and each line contains two space-separated values
      #  matrix = np.array([list(map(float, line.split())) for line in lines])

       # if not matrix.shape[1] == 2:
        #    raise ValueError("Each line in the file should contain two space-separated values.")

       # return matrix


    def construct_mst(self):
        INF = 9999999
        V = len(self.adj_mat)
        selected = [False] * V
        no_edge = 0
        selected[0] = True

        while no_edge < V - 1:
            minimum = INF
            x = 0
            y = 0
            for i in range(V):
                if selected[i]:
                    for j in range(V):
                        if ((not selected[j]) and self.adj_mat[i][j]):  
                            # not in selected and there is an edge
                            if minimum > self.adj_mat[i][j]:
                                minimum = self.adj_mat[i][j]
                                x = i
                                y = j
            print(str(x) + "-" + str(y) + ":" + str(self.adj_mat[x][y]))
            self.mst+=[self.adj_mat[x][y]]
            selected[y] = True
            no_edge += 1
        """
    
        TODO: Given `self.adj_mat`, the adjacency matrix of a connected undirected graph, implement Prim's 
        algorithm to construct an adjacency matrix encoding the minimum spanning tree of `self.adj_mat`. 
            
        `self.adj_mat` is a 2D numpy array of floats. Note that because we assume our input graph is
        undirected, `self.adj_mat` is symmetric. Row i and column j represents the edge weight between
        vertex i and vertex j. An edge weight of zero indicates that no edge exists. 
        
        This function does not return anything. Instead, store the adjacency matrix representation
        of the minimum spanning tree of `self.adj_mat` in `self.mst`. We highly encourage the
        use of priority queues in your implementation. Refer to the heapq module, particularly the 
        `heapify`, `heappop`, and `heappush` functions.

        """
        


#graph_instance=Graph("../data/small.csv")
#Graph.construct_mst(graph_instance)


