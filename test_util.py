'''
Created on Jul 16, 2014

@author: jeromethai
'''

import util
from cvxopt import matrix, spmatrix


def test1():
    """Test bisection"""
    def F(x):
        return 1 + 0.15*(x**4)
    print util.bisection(F, 4.0, 0.0, 5.0)


def test2():
    """Test save_mat"""
    Ms = [matrix(1.0, (3,2)), spmatrix(1.0, range(3), range(3))]
    names = ['M1', 'M2']
    util.save_mat(Ms, names, 'test')


def test3():
    """Read LA_shps_V2"""
    theta = matrix([0.0, 0.0, 0.0, 0.15])
    graph, G, Ids = util.read_shapefile('LA_shps_V2', 'Polynomial', theta, 'LA OSM network')
    print G.number_of_nodes() #return the number of nodes in the graph
    print G.number_of_edges() #return the number of edges b/w two nodes
    print G.size() #return the number of edges
    print G.number_of_selfloops() #return a list of nodes with selfloops
    print graph.numnodes
    print graph.numlinks
    linkIds = [11558, 11425, 20415, 19225, 19529, 5603]
    for id in linkIds: print Ids[id]


def main():
    #test1()
    #test2()
    test3()


if __name__ == '__main__':
    main()