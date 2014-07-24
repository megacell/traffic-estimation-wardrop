'''
Created on Apr 23, 2014

@author: jeromethai
'''
import scipy.linalg as sla
import scipy.io as sio
import numpy as np
from cvxopt import matrix
from numpy.random import normal
import networkx as nx
import Graph as g
import gdal


def place_zeros(M, tol=1e-13):
    """Replace entries in M less than tol by 0.0"""
    for i in range(M.size[0]):
        for j in range(M.size[1]):
            if M[i,j] < tol: M[i,j] = 0.0
    return M


def find_basis(M):
    """Find the indices of the columns of M that form a basis or range(M)"""
    p,l,u = sla.lu(M)
    ind = [i for i in range(u.shape[0]) if u[i,i] != 0.0]
    if u[i,i] == 0:
        for j in range(i+1,u.shape[1]):
            if u[i,j] != 0.0: ind.append(j); break
    return ind


def bisection(F, f, left, right, tol=1e-8):
    """Use bisection to find x such that F(x)=f
    we suppose F is strictly increasing"""
    l, r = left, right
    while r-l>tol:
        if F((l+r)/2.0) < f:
            l = (l+r)/2.0
        else:
            r = (l+r)/2.0
    return (l+r)/2.0


def save_mat(Ms, names, filename):
    """Save matrices in matlab format
    
    Parameters
    __________
    Ms: list of matrices
    names: list of names of the matrices in matlab format
    """
    dict = {names[k] : np.array(matrix(Ms[k])) for k in range(len(Ms))}
    sio.savemat(filename + '.mat', mdict=dict)


def add_noise(A, a, tol=0.1):
    """add gaussian noise to entries of A that are > tol"""
    m,n = A.size
    M = matrix(0.0, (m,n))
    for i in range(m):
        for j in range(n):
            M[i,j] = A[i,j]
            if M[i,j] > tol: M[i,j] = normal(A[i,j], a*A[i,j])
    return M


def create_networkx_graph(graph):
    """Create networkx.DiGraph graph from graph object"""
    G=nx.DiGraph(indedges={})
    G.add_nodes_from(range(1,graph.numnodes))
    G.add_edges_from([(key[0],key[1]) for key in graph.links.keys()])
    i = 0
    for edge in G.edges(): G.graph['indedges'][edge] = i; i+=1
    return G


def read_shapefile(path, delaytype='None', data=None, description=None):
    """Read networkx.DiGraph and return graph.object uses networkx.read_shp
    
    Parameters
    ----------
    path: File, directory, or filename to read by networkx.read_shp
    delaytype: 'None' or 'Polynomial'
    data: if polynomial then data=Theta
    description: description of the graph object
    
    Return value
    ------------
    graph: graph object
    G: networkx.DiGraph object
    IDs: {MATSim link IDs : link ID}
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


if __name__ == '__main__':
    pass
    
    