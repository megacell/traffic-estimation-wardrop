'''
Created on Jul 16, 2014

@author: jeromethai
'''

import util
from cvxopt import matrix, spmatrix
import draw_graph as d
import numpy as np
from multiprocessing import Pool
import time
from cvxopt import matrix, solvers


def test1(x):
    """Test bisection"""
    def F(x):
        return 1 + 0.15*(x**4)
    return util.bisection(F, x, 0.0, 5.0)


def test2():
    """Test save_mat"""
    Ms = [matrix(1.0, (3,2)), spmatrix(1.0, range(3), range(3))]
    names = ['M1', 'M2']
    util.save_mat(Ms, names, 'test')


def test3():
    """Read LA_shps_V2"""
    theta = matrix([0.0, 0.0, 0.0, 0.15])
    graph, G, IDs = util.read_shapefile('LA_shps_V2', 'Polynomial', theta, 'LA OSM network')
    print G.number_of_nodes() #return the number of nodes in the graph
    print G.number_of_edges() #return the number of edges b/w two nodes
    print G.size() #return the number of edges
    print G.number_of_selfloops() #return a list of nodes with selfloops
    print graph.numnodes
    print graph.numlinks
    #d.draw(graph, nodes=False)
    #util.extract_routes_ODs('ODs/grouped_trajectories.csv', 'ODs/grouped_routes.csv', 'ODs/grouped_ODs.csv', IDs, 5)
    util.add_ODs_from_csv(graph, 'ODs/processed_ODs.csv')
    util.add_routes_from_csv(graph, 'ODs/grouped_routes.csv')
    graph.visualize(True)


def test4():
    """Test multiprocessing"""
    start = time.clock()
    pool = Pool(processes=4)
    list = np.linspace(1.0,4.0,2001)
    print pool.map(test1, list)
    t1 = time.clock() - start
    start = time.clock()
    print [test1(x) for x in list]
    t2 = time.clock() - start
    print t1, t2
    
    
def test5(b):
    """Test cvxopt and provide a function using cvxopt package for test6"""
    A = matrix([ [-1.0, -1.0, 0.0, 1.0], [1.0, -1.0, -1.0, -2.0] ])
    c = matrix([ 2.0, 1.0 ])
    sol=solvers.lp(c,A,b)
    return sol['x']
    

def test6():
    """Test cvxopt together with multiprocessing on Python, this works!"""
    b = matrix([ 1.0, -2.0, 0.0, 4.0 ])
    pool = Pool(processes=4)
    print pool.map(test5, [0.5*b, 1.0*b, 1.5*b, 2.0*b])


def main():
    #print test1(4.0)
    #test2()
    #test3()
    #test4()
    #return test5(matrix([ 1.0, -2.0, 0.0, 4.0 ]))
    test6()


if __name__ == '__main__':
    main()