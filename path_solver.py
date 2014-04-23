'''
Created on Apr 22, 2014

@author: jeromethai
'''

from cvxopt import matrix, spmatrix, sparse, solvers
from rank_nullspace import nullspace
from util import place_zeros


def constraints(graph):
    """Construct constraints for feasible path flows
    
    Return value
    ------------
    A: matrix of incidence links-paths
    U: matrix of simplex constraints
    r: matrix of OD flows
    C: matrix -eye(numpaths)
    d: matrix ones(numpaths,1)
    """
    
    I, J = [], []
    for id,path in graph.paths.items():
        for link in path.links: I.append(graph.indlinks[(link.startnode, link.endnode, link.route)]); J.append(graph.indpaths[id])
    A = spmatrix(1.0, I, J)
    
    I, J, r = [], [], matrix(0.0, (graph.numODs,1))
    for id1,od in graph.ODs.items():
        r[graph.indods[id1]] = od.flow
        for id2,path in od.paths.items(): I.append(graph.indods[id1]); J.append(graph.indpaths[id2])
    U = spmatrix(1.0, I, J)
    
    n = graph.numpaths
    C, d = spmatrix(-1.0, range(n), range(n)), matrix(0.0, (n,1))
        
    return A, U, r, C, d


def solver(graph, linkflows, update=True, model='lls', unusedpaths=None, tol=1e-2):
    """Find a feasible path flow
    
    Parameters
    ----------
    graph: graph object
    linkflows: matrix of link flows
    update: if update==True, update path flows in graph
    model: if model=='lls', solve with linear-least-squares
    unusedpaths: ids of unused paths if graph is in UE link flow"""
    
    A, U, r, C, d = constraints(graph)
    
    if model == 'lls': pathflows = solvers.qp(A.trans()*A, -A.trans()*linkflows, C, d, U, r)['x']
        
    if model == 'other': pass
    
    if update == True:
        for id,path in graph.paths.items(): path.flow = pathflows[graph.indpaths[id]]
    
    if unusedpaths is not None:
        valid = True
        for id in unusedpaths:
            if pathflows[graph.indpaths[id]] > tol * graph.ODs[(id[0],id[1])].flow:
                print 'WARNING: path flow ({},{},{}) doesn\'t satisfy UE!'.format(id[0],id[1],id[2]); valid = False; break
            else:
                pathflows[graph.indpaths[id]] = 0.0
        if valid == True: print 'Path flow satisfies UE'
        
    return pathflows


def vec_feas_paths(graph, unusedpaths=None, tol=1e-2):
    """Find a basis of the space of feasible paths"""
    A, U, r, C, d = constraints(graph)
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
 
