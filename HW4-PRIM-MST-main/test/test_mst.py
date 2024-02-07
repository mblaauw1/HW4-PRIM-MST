import numpy as np
import heapq
from typing import Union

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





import pytest
import numpy as np
from mst import Graph
from sklearn.metrics import pairwise_distances


def check_mst(adj_mat: np.ndarray, 
              mst: np.ndarray, 
              expected_weight: int, 
              allowed_error: float = 0.0001):
    """
    
    Helper function to check the correctness of the adjacency matrix encoding an MST.
    Note that because the MST of a graph is not guaranteed to be unique, we cannot 
    simply check for equality against a known MST of a graph. 

    Arguments:
        adj_mat: adjacency matrix of full graph
        mst: adjacency matrix of proposed minimum spanning tree
        expected_weight: weight of the minimum spanning tree of the full graph
        allowed_error: allowed difference between proposed MST weight and `expected_weight`

    TODO: Add additional assertions to ensure the correctness of your MST implementation. For
    example, how many edges should a minimum spanning tree have? Are minimum spanning trees
    always connected? What else can you think of?

    """

    def approx_equal(a, b):
        return abs(a - b) < allowed_error

    total = 0
    print(mst)
    total=sum(mst)
    print(total)
    assert approx_equal(total, expected_weight), 'Proposed MST has incorrect expected weight'


def test_mst_small():
    """
    
    Unit test for the construction of a minimum spanning tree on a small graph.
    
    """
    file_path = '../data/small.csv'
    g = Graph(file_path)
    Graph.construct_mst(g)
    check_mst(g.adj_mat, g.mst, 8)


def test_mst_single_cell_data():
    """
    
    Unit test for the construction of a minimum spanning tree using single cell
    data, taken from the Slingshot R package.

    https://bioconductor.org/packages/release/bioc/html/slingshot.html

    """
    file_path = '../data/slingshot_example.txt'
    coords = np.loadtxt(file_path) # load coordinates of single cells in low-dimensional subspace
    dist_mat = pairwise_distances(coords) # compute pairwise distances to form graph
    g = Graph(dist_mat)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 57.263561605571695)


def test_mst_student(adj_mat: np.ndarray, 
              mst: np.ndarray, 
              allowed_error: float = 0.0001):
    #in a minimal spanning tree for a continuous connected network,
    #the number of edges are fewer than the number of verticies by 
    #edge=n-1
    #where edge=an edge on the network
           #n=node number
    #test will determine if this is true for my minimal spanning trees
    def approx_equal(a, b):
        return abs(a - b) < allowed_error
    edges=len(mst)
    print(edges)
    nodes=len(adj_mat)
    print(nodes)
    assert approx_equal(nodes, edges+1)
    pass

def use_student_test():
    file_path = '../data/slingshot_example.txt'
    coords = np.loadtxt(file_path)
    dist_mat = pairwise_distances(coords) # compute pairwise distances to form graph
    g = Graph(dist_mat)
    g.construct_mst()
    test_mst_student(g.adj_mat, g.mst)




test_mst_small()

test_mst_single_cell_data()

use_student_test()



