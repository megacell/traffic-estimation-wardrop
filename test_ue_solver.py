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
    
    
def los_angeles_ue():
    theta = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
    g1, g2, g3, g4 = los_angeles(theta, 'Polynomial')
    l1 = ue.solver(g1, update=True)
    l2 = ue.solver(g2, update=True)
    l3 = ue.solver(g3, update=True)
    l4 = ue.solver(g4, update=True)
    #d.draw_delays(g1)
    #d.draw_delays(g2)
    #d.draw_delays(g3)
    #d.draw_delays(g4)
    #x1, x2, x3, x4 = np.zeros((g1.numlinks, 1)), np.zeros((g1.numlinks, 1)), np.zeros((g1.numlinks, 1)), np.zeros((g1.numlinks, 1))
    #for i in range(g1.numlinks): x1[i], x2[i], x3[i], x4[i] = l1[i], l2[i], l3[i], l4[i]
    #sio.savemat('linkflows.mat', mdict={'x1': x1, 'x2': x2, 'x3': x3, 'x4': x4})


def main():
    #test1()
    test2()
    #los_angeles_ue()
    

if __name__ == '__main__':
    main()
    