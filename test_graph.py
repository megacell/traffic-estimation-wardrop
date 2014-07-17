'''
Created on Apr 18, 2014

@author: jeromethai
'''

import Graph
import Graph as g
import numpy as np
import scipy.io as sio
from cvxopt import matrix

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


def los_angeles(theta=None, delaytype='None', noisy=False):
    
    data = sio.loadmat('los_angeles_data_2.mat')
        
    if not noisy:
        links = data['links']
        ODs1, ODs2, ODs3, ODs4 = data['ODs1'], data['ODs2'], data['ODs3'], data['ODs4']
    else:
        links = data['links_noisy']
        ODs1, ODs2, ODs3, ODs4 = data['ODs1_noisy'], data['ODs2_noisy'], data['ODs3_noisy'], data['ODs4_noisy']
        
    nodes = data['nodes']
        
    if theta is not None:
        degree = len(theta)
        tmp = links
        links = []
        for startnode, endnode, route, ffdelay, slope in tmp:
            coef = [ffdelay*a*b for a,b in zip(theta, np.power(slope, range(1,degree+1)))]
            links.append((startnode, endnode, route, ffdelay, (ffdelay, slope, coef)))
    
    g1 = g.create_graph_from_list(nodes, links, delaytype, ODs1, 'Map of L.A.')
    g2 = g.create_graph_from_list(nodes, links, delaytype, ODs2, 'Map of L.A.')
    g3 = g.create_graph_from_list(nodes, links, delaytype, ODs3, 'Map of L.A.')
    g4 = g.create_graph_from_list(nodes, links, delaytype, ODs4, 'Map of L.A.')
    return g1, g2, g3, g4


def main():
    #graph = small_example()
    theta = matrix([0.0, 0.0, 0.0, 1.0])
    theta /= np.sum(theta)
    theta *= 0.15
    graph = los_angeles(theta, 'Polynomial')[0]
    graph.visualize(True, True, True, True, True)

if __name__ == '__main__':
    main()
