'''
Created on Apr 18, 2014

@author: jeromethai
'''

import graph as g
import numpy as np
import scipy.io as sio
from cvxopt import matrix
from numpy.random import normal
from util import distance_on_unit_sphere


def small_example():
    
    graph = g.Graph('Small example graph')
    
    graph.add_node((0,1))
    graph.add_node((0,-1))
    graph.add_node((2,0))
    graph.add_node((4,0))
    graph.add_node((6,0))
    
    graph.add_link(1, 3, 1, delayfunc=g.create_delayfunc('Polynomial',(1.0, 1.0, [0.0])))
    graph.add_link(2, 3, 1, delayfunc=g.create_delayfunc('Polynomial',(2.0, 1.0, [0.0])))
    graph.add_link(3, 4, 1, delayfunc=g.create_delayfunc('Polynomial',(2.0, 1.0, [1.0])))
    graph.add_link(3, 4, 2, delayfunc=g.create_delayfunc('Polynomial',(1.0, 1.0, [2.0])))
    graph.add_link(4, 5, 1, delayfunc=g.create_delayfunc('Polynomial',(1.0, 1.0, [0.0])))
    
    graph.add_od(1, 5, 2.0)
    graph.add_od(2, 5, 3.0)
    
    graph.add_path([(1,3,1), (3,4,1), (4,5,1)])
    graph.add_path([(1,3,1), (3,4,2), (4,5,1)])
    graph.add_path([(2,3,1), (3,4,1), (4,5,1)])
    graph.add_path([(2,3,1), (3,4,2), (4,5,1)])
    
    return graph


def los_angeles(parameters=None, delaytype='None', noise=0.0):
    """Generate small map of L.A. with 122 links and 44 modes
    """
    
    data = sio.loadmat('los_angeles_data_2.mat')
        
    links = data['links']
    ODs1, ODs2, ODs3, ODs4 = data['ODs1'], data['ODs2'], data['ODs3'], data['ODs4']
    if noise>0.0:
        ODs1 = [(o, d, normal(f, noise*f)) for o,d,f in ODs1]
        ODs2 = [(o, d, normal(f, noise*f)) for o,d,f in ODs2]
        ODs3 = [(o, d, normal(f, noise*f)) for o,d,f in ODs3]
        ODs4 = [(o, d, normal(f, noise*f)) for o,d,f in ODs4]
        links = [(s, t, r, normal(d, noise*d), c) for s,t,r,d,c in links]
      
    nodes = data['nodes']
    tmp = links
    links = []  
    
    if delaytype=='Polynomial':
        theta = parameters
        degree = len(theta)
        for startnode, endnode, route, ffdelay, slope in tmp:
            coef = [ffdelay*a*b for a,b in zip(theta, np.power(slope, range(1,degree+1)))]
            links.append((startnode, endnode, route, ffdelay, (ffdelay, slope, coef)))
    if delaytype=='Hyperbolic':
        a,b = parameters
        for startnode, endnode, route, ffdelay, slope in tmp:
            k1, k2 = a*ffdelay/slope, b/slope
            links.append((startnode, endnode, route, ffdelay, (ffdelay, slope, k1, k2)))
    if delaytype=='None':
        for startnode, endnode, route, ffdelay, slope in tmp: links.append((startnode, endnode, route, ffdelay, None))
            
    
    g1 = g.create_graph_from_list(nodes, links, delaytype, ODs1, 'Map of L.A.')
    g2 = g.create_graph_from_list(nodes, links, delaytype, ODs2, 'Map of L.A.')
    g3 = g.create_graph_from_list(nodes, links, delaytype, ODs3, 'Map of L.A.')
    g4 = g.create_graph_from_list(nodes, links, delaytype, ODs4, 'Map of L.A.')
    return g1, g2, g3, g4


def los_angeles_2(parameters=None, delaytype='None'):
    """Generate larger map of L.A. with 664 links and 194 nodes
    """
    data = sio.loadmat('los_angeles_data_3.mat')
    nodes, links = data['nodes'], data['links']
    tmp = links
    links = []
    a1 = 60.0*3960.0/75.0
    a2 = 60.0*3960.0/50.0
    
    if delaytype=='None':
        for startnode, endnode, category in tmp:
            arc = distance_on_unit_sphere(nodes[startnode-1][1], nodes[startnode-1][0], nodes[endnode-1][1], nodes[endnode-1][0])
            if category == 1: ffdelay = a1*arc
            if category == 2: ffdelay = a2*arc
            links.append((startnode, endnode, 1, ffdelay, None))
            
    if delaytype=='Polynomial':
        theta = parameters
        degree = len(theta)
        for startnode, endnode, category in tmp:
            arc = distance_on_unit_sphere(nodes[startnode-1][1], nodes[startnode-1][0], nodes[endnode-1][1], nodes[endnode-1][0])
            if category == 1: ffdelay, slope = a1*arc, 5.0
            if category == 2: ffdelay, slope = a2*arc, 1.0
            coef = [ffdelay*a*b for a,b in zip(theta, np.power(slope, range(1,degree+1)))]
            links.append((startnode, endnode, 1, ffdelay, (ffdelay, slope, coef)))

    #print nodes
    #print links
    return g.create_graph_from_list(nodes, links, delaytype, None, 'Larger map of L.A.')


def main():
    #graph = small_example()
    #theta = matrix([0.0, 0.0, 0.0, 0.15])
    #graph = los_angeles(theta, 'Polynomial', 1/15.0)[0]
    #graph.visualize(True, True, True, True, True)
    graph = los_angeles_2(None, 'None')
    graph.visualize(True)

if __name__ == '__main__':
    main()
