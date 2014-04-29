'''
Created on Apr 22, 2014

@author: jeromethai
'''

import numpy as np
import scipy.sparse as sps
import scipy.io as sio
from cvxopt import matrix, spmatrix, sparse, solvers
from rank_nullspace import nullspace
from util import place_zeros


def incidence(graph):
    """Returns matrix of incidence links-paths
    """
    I, J = [], []
    for id,path in graph.paths.items():
        for link in path.links: I.append(graph.indlinks[(link.startnode, link.endnode, link.route)]); J.append(graph.indpaths[id])
    return spmatrix(1.0, I, J)


def simplex(graph):
    """Construct constraints for feasible path flows
    
    Return value
    ------------
    U: matrix of simplex constraints
    r: matrix of OD flows
    """
    I, J, r = [], [], matrix(0.0, (graph.numODs,1))
    for id1,od in graph.ODs.items():
        r[graph.indods[id1]] = od.flow
        for id2,path in od.paths.items(): I.append(graph.indods[id1]); J.append(graph.indpaths[id2])
    return spmatrix(1.0, I, J), r


def solver(graph, linkflows=None, update=True, model='lls', A=None, U=None, r=None, unusedpaths=None, tol=1e-2):
    """Find a feasible path flow
    
    Parameters
    ----------
    graph: graph object
    linkflows: matrix of link flows
    update: if update==True, update path flows in graph
    model: if model=='lls', solve with linear-least-squares
    unusedpaths: ids of unused paths if graph is in UE link flow"""
    
    if A is None: A = incidence(graph)
    if U is None or r is None: U, r = simplex(graph)
    if linkflows is None: print 'Get linkflows from Graph object.'; linkflows = graph.get_linkflows()
        
    n = graph.numpaths
    C, d = spmatrix(-1.0, range(n), range(n)), matrix(0.0, (n,1))
    
    if model == 'lls': pathflows = solvers.qp(A.trans()*A, -A.trans()*linkflows, C, d, U, r)['x']
        
    if model == 'other': pass
    
    if unusedpaths is not None:
        valid = True
        for id in unusedpaths:
            if pathflows[graph.indpaths[id]] > tol * graph.ODs[(id[0],id[1])].flow:
                print 'WARNING: path flow ({},{},{}) doesn\'t satisfy UE!'.format(id[0],id[1],id[2]); valid = False; break
            else:
                pathflows[graph.indpaths[id]] = 0.0
        if valid == True: print 'Path flow satisfies UE'
        
        if update: print 'Update path flows in Graph object.'; graph.update_pathflows(pathflows)
        
    return pathflows


def vec_feas_paths(graph, A=None, U=None, unusedpaths=None, tol=1e-2):
    """Find a basis of the space of feasible paths"""
    if A is None: A = incidence(graph)
    if U is None: U,_ = simplex(graph)
    null = nullspace(matrix(sparse([A, U])))
    ind = range(null.shape[1])
    
    if unusedpaths is not None:
        print 'Trimming unfeasible directions'
        for j in range(null.shape[1]):
            for id in unusedpaths:
                i = graph.indpaths[id]
                if abs(null[i][j]) > tol: ind.remove(j); break
                null[i][j] = 0.0
    
    return place_zeros(null[:,ind])
 

def save_mat(filepath, graph, linkflows):
    """Save sparse matrices for the Path problem in solve_path_data.mat file
    A: incidence matrix link-path
    b: flow on each link
    U: incidence matrix link-OD
    r: OD flows
    """
    I, J = [], []
    b = np.array([linkflows[i] for i in range(len(linkflows))])
    for id,path in graph.paths.items():
        for link in path.links: I.append(graph.indlinks[(link.startnode, link.endnode, link.route)]); J.append(graph.indpaths[id])
    A = sps.coo_matrix(([1]*len(I),(I,J)))
        
    I, J, r = [], [], np.zeros((graph.numODs,1))
    for id1,od in graph.ODs.items():
        r[graph.indods[id1]] = od.flow
        for id2,path in od.paths.items(): I.append(graph.indods[id1]); J.append(graph.indpaths[id2])
    U = sps.coo_matrix(([1]*len(I),(I,J)))
        
    sio.savemat(filepath+'solve_path_data.mat', mdict={'A': A, 'b': b, 'U': U, 'r': r})
    