'''
Created on Apr 21, 2014

@author: jeromethai
'''

import numpy as np
import ue_solver as ue
import draw_graph as d
from test_graph import small_example, los_angeles
from cvxopt import matrix


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
    
    
def los_angeles_ue():
    theta = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
    g4 = los_angeles(theta, 'Polynomial')[3]
    n = g4.numlinks
    g4.add_path([(29,21,1), (21,14,1), (14,34,1), (34,12,1), (12,5,1)])
    g4.add_path([(29,21,1), (21,14,1), (14,13,1), (13,12,1), (12,5,1)])
    g4.add_path([(30,28,1), (28,22,1), (22,15,1), (15,13,1), (13,12,1), (12,5,1)])
    g4.add_path([(30,28,1), (28,23,1), (23,16,1), (16,15,1), (15,13,1), (13,12,1), (12,5,1)])
    l4, x4 = ue.solver(g4, update=True, full=True)
    d.draw_delays(g4, x4[:n])
    d.draw_delays(g4, x4[n:2*n])
    d.draw_delays(g4, x4[2*n:])
    g4.visualize(paths=True)
    

def main():
    #test1()
    #test2()
    los_angeles_ue()
    

if __name__ == '__main__':
    main()
    