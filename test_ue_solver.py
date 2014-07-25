'''
Created on Apr 21, 2014

@author: jeromethai
'''

import numpy as np
import ue_solver as ue
import draw_graph as d
from test_graph import small_example, los_angeles
from cvxopt import matrix, mul


def test1():
    graph = small_example()
    linkflows = ue.solver(graph)
    print 'links\' indices: ', graph.indlinks
    print 'UE flow: '
    print linkflows


def test2():
    theta = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
    graph = los_angeles(theta, 'Polynomial')[0]
    C,ind = ue.get_nodelink_incidence(graph)
    print C
    d = ue.get_demands(graph, ind, 22)
    print d
    Aeq, beq = ue.constraints(graph)
    print Aeq.size
    print beq.size
    
    
def test3():
    theta = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
    g4 = los_angeles(theta, 'Polynomial')[3]
    n = g4.numlinks
    g4.add_path_from_nodes([29,21,14,34,12,5])
    g4.add_path_from_nodes([29,21,14,13,12,5])
    g4.add_path_from_nodes([30,28,22,15,13,12,5])
    g4.add_path_from_nodes([30,28,23,16,15,13,12,5])
    l4, x4 = ue.solver(g4, update=True, full=True)
    #d.draw_delays(g4, x4[:n])
    #d.draw_delays(g4, x4[n:2*n])
    #d.draw_delays(g4, x4[2*n:])
    print l4
    #g4.visualize(paths=True)
    

def main():
    #test1()
    #test2()
    test3()
    

if __name__ == '__main__':
    main()
    