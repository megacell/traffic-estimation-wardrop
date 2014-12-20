# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 12:10:30 2014

@author: hugo
"""
from cvxopt import matrix
import read_csv as r
import numpy as np
import graph as g
from util import closest_node

def osm_network(parameters=None, delaytype='None', box=None):
    
    latmin, latmax, lngmin, lngmax = box[0], box[1], box[2], box[3]  
    nodes, links = r.create_nodes_links_from_box(latmin, latmax, lngmin, lngmax)
    links = r.kill_double_links(links)   
    '''        
    for link in links:#To simulate an incident
        if (link[0]==36 and link[1]==38):
            link[4]=10*link[4]
    '''  
    tmp = links
    links=[]
    if delaytype=='Polynomial':
        theta = parameters
        degree = len(theta)
        for startnode, endnode, route, ffdelay, slope in tmp:
            coef = [ffdelay*a*b for a,b in zip(theta, np.power(slope, range(1,degree+1)))]
            links.append((startnode, endnode, route, ffdelay, (ffdelay, slope, coef)))      
            
    dest1 = 205
    dest2 = 207
    ODs=[]
    ODs+=create_linear_ODs(34.176831, -118.131931, 34.157088, -118.018977, 10, dest1, nodes, 3000.0)
    ODs+=create_linear_ODs(34.132226, -118.10755, 34.128247, -118.024299, 10, dest1, nodes, 3000.0)
    #ODs+=create_linear_ODs(34.176831, -118.131931, 34.157088, -118.018977, 10, dest2, nodes, 3000.0)
    #ODs+=create_linear_ODs(34.132226, -118.10755, 34.128247, -118.024299, 10, dest2, nodes, 3000.0)
    #print links[0:50]
    #ODs=[[1, 2 ,5.0]]
    print ODs
    g1 = g.create_graph_from_list(nodes, links, delaytype, ODs, 'Map of L.A.')
    return g1
    
def main():
    theta = matrix([0.0, 0.0, 0.0, 0.15])
    #theta = matrix([0.0, 0.15])
    #graph = xml_network(theta, 'Polynomial') # I got rid of the 0.15 .... Hugo
    graph.visualize(True, True, True, True, True)

if __name__ == '__main__':
    main()
    
def create_linear_ODs(lat1, lng1, lat2, lng2, nb, dest, nodes, value_tot):
    ODs = []
    for i in range(nb):
        lat, lng = lat1 + i/float(nb)*(lat2-lat1), lng1 + i/float(nb)*(lng2-lng1)
        ODs.append([closest_node(lat, lng, nodes),dest,float(value_tot)/nb])
    return ODs