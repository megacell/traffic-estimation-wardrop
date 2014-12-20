# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 18:20:27 2014

@author: hugo
"""
import numpy as np
from cvxopt import matrix
import graph as g
from util import distance_on_unit_sphere
   
def Get_nodes_from_csv():
    data = np.genfromtxt('Data/nodescsv.csv', delimiter = ',', skiprows=1)
    Nodes_table=[]
    for i in range(len(data)):
        node=[data[i][0], data[i][1], data[i][2]]
        Nodes_table.append(node)
    return Nodes_table

            
def Get_links_from_csv(Nodes):
    data = np.genfromtxt('Data/links_qgis.csv', delimiter = ',', skiprows=1)
    Links_table=[]
    theta = matrix([0.0, 0.0, 0.0, 0.15])
    degree=len(theta)
    for i in range(len(data)):
        startnode, endnode, type= int(data[i][0]), int(data[i][1]), data[i][2]
        length = distance_on_unit_sphere(Nodes[startnode-1][2], Nodes[startnode-1][1], Nodes[endnode-1][2], Nodes[endnode-1][1])
        if type==1:
            ffdelay=length / 32
            slope = 0.2
        else: ffdelay, slope = length / 17, 0.75
        coef = [ffdelay*a*b for a,b in zip(theta, np.power(slope, range(1,degree+1)))]
        Links_table.append((startnode, endnode, 1, ffdelay, (ffdelay, slope, coef))) 
    return Links_table

def tests():
    Nodes = Get_nodes_from_csv()
    Links = Get_links_from_csv(Nodes)
    ODs = [[68,70,5.0]]
    g1 = g.create_graph_from_list(Nodes, Links, 'Polynomial', ODs, 'Map of L.A.')
    return g1    
              
def main():
    #tests()
    return 0
    
if __name__ == '__main__':
    main()