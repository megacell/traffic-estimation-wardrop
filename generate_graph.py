'''
Created on Apr 18, 2014

@author: jeromethai
'''
import Graph as g
import numpy as np
import scipy.io as sio
from cvxopt import matrix
from numpy.random import normal
from util import distance_on_unit_sphere
import draw_graph as d
from get_ODs_from_csv import Create_ODs_nodes_unique

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


def los_angeles(parameters=None, delaytype='None', noise=0.0, path=None):
    """Generate small map of L.A. with 122 links and 44 modes
    """

    if not path:
        path = 'los_angeles_data_2.mat'
    data = sio.loadmat(path)

        
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
    nodes = np.genfromtxt('LA_medium_data/nodes_LA_toy.csv', delimiter = ',', skiprows = 1)
    nodes = nodes[:,1:3]
    links = np.genfromtxt('LA_medium_data/links_qgis_cap.csv', delimiter = ',', skiprows = 1)
    tmp = links
    links = []
    speed_limit_freeway = 33.33 #unit: m/s
    dict_cap2speed ={600:12.5, 1000:16.67, 2000:16.67, 4000:16.67, 5000:16.67, 1500:16.67, 3000:22.22, 6000:22.22, 9000:22.22, 4500:22.22, 7500:22.22, 10500:22.22}
    
    if delaytype=='None':
        for startnode, endnode, category in tmp:
            arc = distance_on_unit_sphere(nodes[startnode-1][2], nodes[startnode-1][1], nodes[endnode-1][2], nodes[endnode-1][1])
            if category == 1: ffdelay = arc/speed_limit_freeway
            if category == 2: ffdelay = arc/16.67
            if category !=0: links.append((startnode, endnode, 1, ffdelay, None))
            
    if delaytype=='Polynomial':
        theta = parameters
        degree = len(theta)
        for startnode, endnode, category, cap in tmp:
            arc = distance_on_unit_sphere(nodes[startnode-1][1], nodes[startnode-1][0], nodes[endnode-1][1], nodes[endnode-1][0])
            if category == 1: ffdelay, slope = arc/speed_limit_freeway, 1/cap
            if category == 2: 
                if dict_cap2speed.has_key(cap):           
                    ffdelay, slope = arc/dict_cap2speed[cap], 1/cap
                else : ffdelay, slope = arc/16.67, 1/cap
            coef = [ffdelay*a*b for a,b in zip(theta, np.power(slope, range(1,degree+1)))]
            links.append((startnode, endnode, 1, ffdelay, (ffdelay, slope, coef)))
    
    dest1 = 50
    dest2 = 100
    ODs=[]
    
    #ODs+=create_linear_ODs(34.044801, -117.831116, 33.955998, -118.309021, 10, 10, nodes, 6000.0)
    #ODs+=create_linear_ODs(34.162493, -118.301468, 34.106226, -117.903214, 2, 20, nodes, 6000.0)    
    #ODs+=create_linear_ODs(34.044801, -117.831116, 33.955998, -118.309021, 2, 30, nodes, 6000.0)
    #ODs+=create_linear_ODs(34.162493, -118.301468, 34.106226, -117.903214, 2, 40, nodes, 6000.0)    
    #ODs+=create_linear_ODs(34.044801, -117.831116, 33.955998, -118.309021, 2, 50, nodes, 6000.0)
    #ODs+=create_linear_ODs(34.162493, -118.301468, 34.106226, -117.903214, 2, 60, nodes, 6000.0)    
    #ODs+=create_linear_ODs(34.044801, -117.831116, 33.955998, -118.309021, 2, 70, nodes, 6000.0)
    #ODs+=create_linear_ODs(34.162493, -118.301468, 34.106226, -117.903214, 2, 80, nodes, 6000.0)        
        
    #ODs = [[  57. ,   50.  ,   4. ],[  20. ,  162.   , 56.4]]
    ODs = Create_ODs_nodes_unique(nodes)
    ODs = ODs[1:5]
    print ODs
    return g.create_graph_from_list(nodes, links, delaytype, ODs, 'Larger map of L.A.')




def main():
    theta = matrix([0.0, 0.0, 0.0, 0.15])
    #graph = los_angeles_2(theta, 'Polynomial', 1/15.0)[0]
    graph = los_angeles_2(theta, 'Polynomial')
    #graph.visualize(True, True, True, True, True)
    d.draw(graph, nodes=False)

if __name__ == '__main__':
    main()
