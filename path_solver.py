'''
Created on Apr 22, 2014

@author: jeromethai
'''

from cvxopt import matrix, spmatrix


def index_paths(graph):
    """Index paths of the graph with integers from 0 to numpaths-1"""
    return {graph.paths.keys()[i]: i for i in range(graph.numpaths)}


def index_ods(graph):
    """Index ODs of the graph with integers from 0 to numODs-1"""
    return {graph.ODs.keys()[i]: i for i in range(graph.numODs)}


def constraints(graph, linkflows, indlinks):
    """Construct constraints for feasible path flows
    
    Return value
    ------------
    indpaths: dict of path id to index
    indods: dict of od id to index
    A: matrix of incidence links-paths
    b: matrix of link flows
     U: matrix of simplex constraints
    r: matrix of OD flows
    C: matrix -eye(numpaths)
    d: matrix ones(numpaths,1)
    """
    
    indpaths = index_paths(graph)
    I, J = [], []
    for id,path in graph.paths.items():
        for link in path.links: I.append(indlinks[(link.startnode, link.endnode, link.route)]); J.append(indpaths[id])
    A = spmatrix(1.0, I, J)
    b = linkflows
    
    indods = index_ods(graph)
    I, J, r = [], [], matrix(0.0, (graph.numODs,1))
    for id1,od in graph.ODs.items():
        r[indods[id1]] = od.flow
        for id2,path in od.paths.items(): I.append(indods[id1]); J.append(indpaths[id2])
    U = spmatrix(1.0, I, J)
    
    n = graph.numpaths
    C, d = spmatrix(-1.0, range(n), range(n)), matrix(0.0, (n,1))
        
    return indpaths, indods, A, b, U, r, C, d


def solver(graph, linkflows, indlinks, update=False, model='lls'):
    """Find a feasible path flow
    if update==True: update path flows in graph"""
    
    indpaths, indods, A, b, U, r, C, d = constraints(graph, linkflows, indlinks)
    
    if model == 'lls':
        from cvxopt.solvers import qp
        pathflows = qp(A.trans()*A, -A.trans()*b, C, d, U, r)['x']
        
    if model == 'other':
        pass
    
    return indpaths, pathflows