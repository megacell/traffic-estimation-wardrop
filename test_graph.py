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
    
    graph.add_link(1, 3, 1, delayfunc=g.create_delayfunc('Affine',(1.0, 0.0)))
    graph.add_link(2, 3, 1, delayfunc=g.create_delayfunc('Affine',(2.0, 0.0)))
    graph.add_link(3, 4, 1, delayfunc=g.create_delayfunc('Affine',(2.0, 1.0)))
    graph.add_link(3, 4, 2, delayfunc=g.create_delayfunc('Affine',(1.0, 2.0)))
    graph.add_link(4, 5, 1, delayfunc=g.create_delayfunc('Affine',(1.0, 0.0)))
    
    graph.add_od(1, 5, 2.0)
    graph.add_od(2, 5, 3.0)
    
    graph.add_path([(1,3,1), (3,4,1), (4,5,1)])
    graph.add_path([(1,3,1), (3,4,2), (4,5,1)])
    graph.add_path([(2,3,1), (3,4,1), (4,5,1)])
    graph.add_path([(2,3,1), (3,4,2), (4,5,1)])
    
    return graph


def small_grid(od_flows, delaytype='Affine', theta=None):
    """Creates a small grid with fixed geometry, fixed affine/polynomial latency functions, and fixed OD pairs
    variable OD flows
    
    Parameters
    ----------
    od_flows: OD demands
    delaytype: type  of the delay functions
    theta: if delaytype = 'Polynomial', 
            then delay at link i is D_i(x) = ffdelays[i] + sum^{degree}_{k=1} theta[k-1]*(slopes[i]*x)^k
    """
    
    ffdelays = [1.0, 1.0, 1.0, 3.0, 2.0, 2.0, 3.0, 1.0]
    #slopes = [4.0, 4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 4.0]
    slopes = [1.0, 1.0, 0.25, 0.25, 0.25, 0.25, 0.25, 1.0]
    
    if delaytype == 'Affine': data = zip(ffdelays, slopes)
    elif delaytype == 'Polynomial':
        degree = len(theta)
        coefs = []
        for j in range(8):
            coefs.append([ffdelays[j]*a*b for a,b in zip(theta, np.power(slopes[j], range(1,degree+1)))])
        data = zip(ffdelays, slopes, coefs)
    
    grid = g.create_grid(2, 3, outdown=[1,1,1],
        outddelay=[[data[0]], [data[1]], [data[2]]], 
        inright=[1,1,0,1,1],
        inrdelay=[[data[3]], [data[4]], [(0.0, 0.0)], [data[5]], [data[6]]],
        delaytype=delaytype)
    
    
    grid.add_link(5, 2, 1, delayfunc=g.create_delayfunc(delaytype,data[7]))
        
    grid.add_od(3, 6, od_flows[0])
    grid.add_od(3, 4, od_flows[1])
    grid.add_od(3, 2, od_flows[2])
    grid.add_od(2, 4, od_flows[3])
    
    #grid.add_od(10,3) # error
    #grid.add_od(0,3) # error
    #grid.add_od(2,4) # error
    
    # paths for OD (3,6)
    grid.add_path([(3,6,1)])
    
    #paths for OD (3,4)
    grid.add_path([(3,2,1), (2,1,1), (1,4,1)])
    grid.add_path([(3,2,1), (2,5,1), (5,4,1)])
    grid.add_path([(3,6,1), (6,5,1), (5,4,1)])
    grid.add_path([(3,6,1), (6,5,1), (5,2,1), (2,1,1), (1,4,1)])
    
    #grid.add_path([(3,2,1), (2,5,1)]) # error
    #grid.add_path([(3,2,1), (2,5,1), (5,4,1)]) # error
    
    # paths for OD (3,2)
    grid.add_path([(3,2,1)])
    grid.add_path([(3,6,1), (6,5,1), (5,2,1)])
    
    # paths for OD (2,4)
    grid.add_path([(2,1,1), (1,4,1)])
    grid.add_path([(2,5,1), (5,4,1)])
    return grid


def los_angeles(theta=None, delaytype='None', noisy=False):
    
    data = sio.loadmat('los_angeles_data.mat')
        
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
    #graph = small_grid([3.0, 3.0, 1.0, 1.0])
    #graph = small_grid([3.0, 3.0, 1.0, 1.0], 'Polynomial', [1.0, 2.0, 3.0])
    #graph.visualize(True, True, True, True, True)
    theta = matrix([0.0, 0.0, 0.0, 1.0])
    theta /= np.sum(theta)
    theta *= 0.15
    graph = los_angeles(theta, 'Polynomial')
    graph.visualize(True, True, True, True, True)
    

if __name__ == '__main__':
    main()
    
    
