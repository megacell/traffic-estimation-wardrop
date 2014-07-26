'''
Created on Jul 23, 2014

@author: jeromethai
'''

import ue_solver as ue
from cvxopt import matrix, spmatrix, spdiag, solvers, div

def linkpath_incidence(graph):
    """Returns matrix of incidence link-path
    """
    I, J, entries = [], [], []
    for id1,link in graph.links.items():
        if len(link.paths) > 0:
            for id2 in link.paths.keys():
                I.append(graph.indlinks[id1])
                J.append(graph.indpaths[id2])
                entries.append(1.0)
        else:
            I.append(graph.indlinks[id1])
            J.append(0)
            entries.append(0.0)
    return spmatrix(entries, I, J)


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
        for id2,path in od.paths.items():
            I.append(graph.indods[id1])
            J.append(graph.indpaths[id2])
    return spmatrix(1.0, I, J), r


def solver_init(U,r):
    n,m = U.size
    x0 = matrix(0.0, (m,1))
    for i in range(n):
        k = sum(U[i,:])
        if k > 0.0:
            for j in range(m):
                if U[i,j] > 0.0: x0[j] = r[i]/k
    return x0


def solver(graph, update=False, data=None):
    """Solve for the UE equilibrium using link-path formulation
    
    Parameters
    ----------
    graph: graph object
    update: if update==True, update path flows in graph
    data: (P,U,r) 
            P: link-path incidence matrix
            U,r: simplex constraints 
    """
    type = graph.links.values()[0].delayfunc.type
    if data is None:
        P = linkpath_incidence(graph)
        U,r = simplex(graph)
    else: (P,U,r) = data
    m = graph.numpaths
    A, b = spmatrix(-1.0, range(m), range(m)), matrix(0.0, (m,1))
    ffdelays = graph.get_ffdelays()
    if type == 'Polynomial':
        coefs = graph.get_coefs()
        coefs_i = coefs * spdiag([1.0/(j+2) for j in range(coefs.size[1])])
        parameters = matrix([[ffdelays], [coefs_i]])
        G = ue.objective_poly
    if type == 'Hyperbolic':
        ks = graph.get_ks()
        parameters = matrix([[ffdelays-div(ks[:,0],ks[:,1])], [ks]])
        G = ue.objective_hyper
    def F(x=None, z=None):
        if x is None: return 0, solver_init(U,r)
        if z is None:
            f, Df = G(P*x, z, parameters, 1)
            return f, Df*P
        f, Df, H = G(P*x, z, parameters, 1)
        return f, Df*P, P.T*H*P    
    x = solvers.cp(F, G=A, h=b, A=U, b=r)['x']
    if update:
        print 'Update link flows, delays in Graph.'; graph.update_linkflows_linkdelays(A*x)
        print 'Update path delays in Graph.'; graph.update_pathdelays()
        print 'Update path flows in Graph object.'; graph.update_pathflows(x)
    return x
