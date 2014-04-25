'''
Created on Apr 23, 2014

@author: jeromethai
'''

import numpy as np
import path_solver as path
from cvxopt import matrix


def remove_meas(graph, link_ids, A=None, linkflows=None):
    """Remove measurements from links with ids link_ids
    
    Parameters
    ----------
    graph: graph object
    link_ids: list of ids [id1 , id2, id3, ...] with ids of the form id = (link.startnode, link.endnode, link.route)
    A: matrix of incidence links-paths
    linkflows: matrix of link flows
    """
    if A is None: A = path.incidence(graph)
    if linkflows is None: print 'Get linkflows from Graph object.'; linkflows = graph.get_linkflows()

    ind = range(graph.numlinks)
    for linkid in link_ids: ind.remove(graph.indlinks[linkid])
    return ind, A[ind,:], linkflows[ind,:]


def remove_meas_rand(graph, k, A=None, linkflows=None):
    """Remove randomly k measurements from links"""
    if A is None: A = path.incidence(graph)
    if linkflows is None: print 'Get linkflows from Graph object.'; linkflows = graph.get_linkflows()
        
    n = graph.numlinks; ind = range(n)
    tmp = np.random.permutation(n)[range(k)]
    for i in tmp: ind.remove(i)
        
    return ind, A[ind,:], linkflows[ind,:]


def error(graph, link_ids=None, k=None, A=None, linkflows=None):
    """Compute error in link flows when measurements are missing"""
    if A is None: A = path.incidence(graph)
    if linkflows is None: print 'Get linkflows from Graph object.'; linkflows = graph.get_linkflows()

    return