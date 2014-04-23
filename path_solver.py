'''
Created on Apr 22, 2014

@author: jeromethai
'''

<<<<<<< HEAD
from cvxopt import matrix, spmatrix


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


def solver(graph, linkflows, update=False, model='lls'):
    """Find a feasible path flow
    
    Parameters
    ----------
    graph: graph object
    linkflows: matrix of link flows
    update: if update==True, update path flows in graph
    model: if model=='lls', solve with linear-least-squares """
    
    A, U, r, C, d = constraints(graph)
    
    if model == 'lls':
        from cvxopt.solvers import qp
        pathflows = qp(A.trans()*A, -A.trans()*linkflows, C, d, U, r)['x']
        
    if model == 'other':
        pass
    
    return pathflows
=======
def constraints(graph, linkflows, indlinks):
    """Construct constraints for feasible path flows"""
    return indpaths, C, d, A, d, Aeq, beq


def solver(graph, linkflows, indlinks, update=False, type='lls'):
    """Find a feasible path flow"""
    return indpaths, pathflows
>>>>>>> parent of 47bc38d... complete path_solver.py and test_path_solver
