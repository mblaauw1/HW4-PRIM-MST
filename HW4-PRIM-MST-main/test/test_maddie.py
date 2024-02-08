import numpy as np
from typing import Union
import pytest
#from mst import Graph
from sklearn.metrics import pairwise_distances


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
    file_path ='../HW4-PRIM-MST-main/data/small.csv'
    g = Graph(file_path)
    Graph.construct_mst(g)
    check_mst(g.adj_mat, g.mst, 8)


def test_mst_single_cell_data():
    """
    
    Unit test for the construction of a minimum spanning tree using single cell
    data, taken from the Slingshot R package.

    https://bioconductor.org/packages/release/bioc/html/slingshot.html

    """
    file_path = '../HW4-PRIM-MST-main/data/slingshot_example.txt'
    coords = np.loadtxt(file_path) # load coordinates of single cells in low-dimensional subspace
    dist_mat = pairwise_distances(coords) # compute pairwise distances to form graph
    g = Graph(dist_mat)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 57.263561605571695)



def check_verticies_vs_conections(adj_mat: np.ndarray, 
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

def test_mst_student():
    file_path = '../HW4-PRIM-MST-main/data/slingshot_example.txt'
    coords = np.loadtxt(file_path)
    dist_mat = pairwise_distances(coords) # compute pairwise distances to form graph
    g = Graph(dist_mat)
    g.construct_mst()
    check_verticies_vs_conections(g.adj_mat, g.mst)




#test_mst_small()

#test_mst_single_cell_data()

#use_student_test()



