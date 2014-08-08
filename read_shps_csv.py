'''
Created on Aug 7, 2014

@author: jeromethai
'''

import networkx as nx
import Graph as g
import gdal
import csv
from cvxopt import matrix
import numpy as np


def read_shapefile(path, delaytype='None', data=None, description=None):
    """Read networkx.DiGraph and return graph.object uses networkx.read_shp
    
    Parameters
    ----------
    path: File, directory, or filename to be read by networkx.read_shp
    delaytype: 'None' or 'Polynomial'
    data: if polynomial then data=Theta
    description: description of the graph object
    
    Return value
    ------------
    graph: graph object
    G: networkx.DiGraph object
    IDs: {MATSim link ID : graph link ID}
    """
    G = nx.read_shp(path)
    nodes, edges = G.nodes(), G.edges(data=True)
    d = {key:i+1 for i,key in enumerate(nodes)}
    IDs = {int(e[2]['ID']): (d[e[0]], d[e[1]],1) for e in edges}
    if delaytype == 'None':
        links = [(d[e[0]], d[e[1]], 1, (e[2]['length']/e[2]['freespeed'])/60.0, None) for e in edges]
    if delaytype == 'Polynomial':
        links, degree = [], len(data)
        for e in edges:
            ffdelay = (e[2]['length'] / e[2]['freespeed'])/60.0 # in minutes
            slope = 1/(e[2]['capacity']/2000.0)
            coef = [ffdelay*a*b for a,b in zip(data, np.power(slope, range(1,degree+1)))]
            links.append((d[e[0]], d[e[1]], 1, ffdelay, (ffdelay, slope, coef)))
    graph = g.create_graph_from_list(nodes, links, delaytype, description=description)
    return graph, G, IDs
                

def extract_routes_ODs(pathin, pathout1, pathout2, IDs, k):
    """Extract routes from .csv file
    
    Parameters
    ----------
    pathin: path to the input file
    pathout1: path to the output file
    IDs: {MATSim link ID: graph link ID}
    k: minimum number of routes in the group
    """
    file = open(pathout1, "w")
    with open(pathin, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            a, b, c = int(row[1][1:-1].partition(', ')[0]), int(row[1][1:-1].rpartition(', ')[-1]), int(row[2])
            if a!=b and c>=k:
                try:
                    iter, node_ids = 0, []
                    for i in row[1][1:-1].split(', '):
                        id = int(i)
                        if iter==0: node_ids.append(IDs[id][0]); node_ids.append(IDs[id][1])
                        if iter>0: node_ids.append(IDs[id][1])
                        iter += 1
                    if node_ids[0] != node_ids[-1]:
                        for id in node_ids: file.write('{},'.format(id))
                        file.write('{}\n'.format(c))
                except KeyError, e:
                    print 'Oops!  That was no valid MATSim id. KeyError: {}'.format(e)
    
    file = open(pathout2, "w")
    with open(pathout1, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            file.write('{},{},{}\n'.format(int(row[0]),int(row[-2]),int(row[-1])))
    

def add_ODs_from_csv(graph, pathin):
    """Add ODs to graph object"""
    with open(pathin, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            graph.add_od(int(row[0]),int(row[1]),float(row[2]))


def add_routes_from_csv(graph, pathin):
    """Add routes to graph object"""
    with open(pathin, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            node_ids = [int(id) for id in row]
            node_ids = node_ids[:-1]
            graph.add_path_from_nodes(node_ids)


def main():
    """Read LA_shps_V2"""
    theta = matrix([0.0, 0.0, 0.0, 0.15])
    graph, G, IDs = read_shapefile('LA_shps_V2', 'Polynomial', theta, 'LA OSM network')
    print G.number_of_nodes() #return the number of nodes in the graph
    print G.number_of_edges() #return the number of edges b/w two nodes
    print G.size() #return the number of edges
    print G.number_of_selfloops() #return a list of nodes with selfloops
    print graph.numnodes
    print graph.numlinks
    #d.draw(graph, nodes=False)
    #util.extract_routes_ODs('ODs/grouped_trajectories.csv', 'ODs/grouped_routes.csv', 'ODs/grouped_ODs.csv', IDs, 5)
    add_ODs_from_csv(graph, 'ODs/processed_ODs.csv')
    add_routes_from_csv(graph, 'ODs/grouped_routes.csv')
    graph.visualize(True)


if __name__ == '__main__':
    main()