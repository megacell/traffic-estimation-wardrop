'''
Created on Apr 21, 2014

@author: jeromethai
'''

import numpy as np
import ue_solver as ue
import draw_graph as d
from generate_graph import small_example, los_angeles, los_angeles_2
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
    #g = los_angeles(theta, delaytype)[3]
    g = los_angeles_2(theta, delaytype)
    n = g.numlinks
    l, x = ue.solver(g, update=True, full=True)
    d.draw_delays(g, x[:n])
    d.draw_delays(g, x[n:2*n])
    d.draw_delays(g, x[2*n:])
    #print l
    print max(mul(l,g.get_slopes()))
    print 'cost UE:', sum([link.delay*link.flow for link in g.links.values()])
    l2, x2 = ue.solver(g, update=True, full=True, SO=True)
    #d.draw_delays(g, x2[:n])
    #d.draw_delays(g, x2[n:2*n])
    #d.draw_delays(g, x2[2*n:])
    #print l2
    print max(mul(l2,g.get_slopes()))
    print 'cost SO:', sum([link.delay*link.flow for link in g.links.values()])


def main():
    #test1()
    test2('Polynomial')
    #test2('Hyperbolic')


if __name__ == '__main__':
    main()
    