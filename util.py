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
    

if __name__ == '__main__':
    pass
    
    