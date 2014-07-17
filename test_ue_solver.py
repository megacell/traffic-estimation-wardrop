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
    g1, g2, g3, g4 = los_angeles(theta, 'Polynomial')
    n = g1.numlinks
    l1, x1 = ue.solver(g1, update=True, full=True)
    l2, x2 = ue.solver(g2, update=True, full=True)
    l3, x3 = ue.solver(g3, update=True, full=True)
    l4, x4 = ue.solver(g4, update=True, full=True)
    d.draw_delays(g1, x1[:n])
    d.draw_delays(g1, x1[n:])
    d.draw_delays(g2, x2[:n])
    d.draw_delays(g2, x2[n:])
    d.draw_delays(g3, x3[:n])
    d.draw_delays(g3, x3[n:])
    d.draw_delays(g4, x4[:n])
    d.draw_delays(g4, x4[n:])
    #print l4[g1.indlinks[(12,5,1)]]
    #print l4[g1.indlinks[(6,5,1)]]


def main():
    #test1()
    #test2()
    los_angeles_ue()
    

if __name__ == '__main__':
    main()
    