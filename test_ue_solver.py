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

    
def los_angeles_ue():
    theta = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
    g1, g2, g3, g4 = los_angeles(theta, 'Polynomial')
    l1 = ue.solver(g1, update=True)
    l2 = ue.solver(g2, update=True)
    l3 = ue.solver(g3, update=True)
    l4 = ue.solver(g4, update=True)
    d.draw_delays(g1)
    d.draw_delays(g2)
    d.draw_delays(g3)
    d.draw_delays(g4)
    

def main():
    #test1()
    los_angeles_ue()
    

if __name__ == '__main__':
    main()
    