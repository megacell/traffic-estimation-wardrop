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


def test2(delaytype):
    if delaytype == 'Polynomial': theta = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
    if delaytype == 'Hyperbolic': theta = (3.5, 3.0)
    g = los_angeles(theta, delaytype)[3]
    n = g.numlinks
    g.add_path_from_nodes([29,21,14,34,12,5])
    g.add_path_from_nodes([29,21,14,13,12,5])
    g.add_path_from_nodes([30,28,22,15,13,12,5])
    g.add_path_from_nodes([30,28,23,16,15,13,12,5])
    l, x = ue.solver(g, update=True, full=True)
    #d.draw_delays(g, x[:n])
    #d.draw_delays(g, x[n:2*n])
    #d.draw_delays(g, x[2*n:])
    print l
    print max(mul(l,g.get_slopes()))
    #g.visualize(paths=True)


def main():
    #test1()
    test2('Polynomial')
    #test2('Hyperbolic')
    

if __name__ == '__main__':
    main()
    