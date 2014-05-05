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


def error_linkflows(graph, link_ids=None, k=None, A=None, linkflows=None, U=None, r=None, ord=None, relative=True):
    """Compute (relative) error in link flows when measurements are missing
    
    Parameters
    ----------
    graph: graph object
    link_ids: list of ids [id1 , id2, id3, ...] with ids of the form id = (link.startnode, link.endnode, link.route)
    k: number of measurements removed randomly
    A: matrix of incidence links-paths
    linkflows: matrix of link flows
    U: matrix of simplex constraints
    r: matrix of OD flows
    ord: order of the norm in the error
    """
    if (link_ids is None) ^ (k is not None): print 'ERROR: must have exactly one of link_ids, k arguments.'; return
    if A is None: A = path.incidence(graph)
    if linkflows is None: print 'Get linkflows from Graph object.'; linkflows = graph.get_linkflows()
    if U is None or r is None: U, r = path.simplex(graph)
    
    if link_ids is not None: ind, misA, misflows = remove_meas(graph, link_ids, A, linkflows)
    else: ind, misA, misflows = remove_meas_rand(graph, k, A, linkflows)
    
    linkflows2 = A*path.solver(graph, misflows, False, 'lls', misA, U, r)
    error = np.linalg.norm(linkflows-linkflows2, ord)
    if relative: error /= np.linalg.norm(linkflows, ord)
    return error, linkflows2


def avg_error(graph, ks, trials, A=None, linkflows=None, U=None, r=None, ord=None, relative=True):
    """Compute average errors in link flows when k measurements are randomly missing
    
    Parameters
    ----------
    graph: graph object
    ks: list of number of measurements removed randomly
    A: matrix of incidence links-paths
    linkflows: matrix of link flows
    U: matrix of simplex constraints
    r: matrix of OD flows
    ord: order of the norm in the error
    """
    if A is None: A = path.incidence(graph)
    if linkflows is None: print 'Get linkflows from Graph object.'; linkflows = graph.get_linkflows()
    if U is None or r is None: U, r = path.simplex(graph)
    avg_errors = []
    for i in range(len(ks)):
        avg_error = 0
        for j in range(trials): avg_error += error_linkflows(graph, None, ks[i], A, linkflows, U, r, ord, relative)[0]
        avg_errors.append(avg_error/trials)
    return avg_errors
    
    