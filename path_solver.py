'''
Created on Jul 23, 2014

@author: jeromethai
'''

import ue_solver as ue
from cvxopt import matrix, spmatrix, spdiag, solvers, div
import numpy.random as ra


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


def solver_init(U,r, random=False):
    """Initialize with a feasible point
    
    Parameters:
    ---------
    U: matrix of simplex constraints
    r: matrix of OD flows
    random: if true
    """
    n,m = U.size
    x0 = matrix(0.0, (m,1))
    for i in range(n):
        k = int(sum(U[i,:]))
        if random:
            tmp = ra.rand(k)
            tmp = r[i]*(tmp/sum(tmp))
        else: tmp = [r[i]/k]*k
        k = 0
        for j in range(m):
            if U[i,j] > 0.0: x0[j] = tmp[k]; k += 1
    return x0


def solver(graph, update=False, data=None, SO=False, random=False):
    """Solve for the UE equilibrium using link-path formulation
    
    Parameters
    ----------
    graph: graph object
    update: if True, update link and path flows in graph
    data: (P,U,r) 
            P: link-path incidence matrix
            U,r: simplex constraints 
    SO: if True compute SO
    random: if True, initialize with a random feasible point
    """
    type = graph.links.values()[0].delayfunc.type
    if data is None:
        P = linkpath_incidence(graph)
        U,r = simplex(graph)
    else: P,U,r = data
    m = graph.numpaths
    A, b = spmatrix(-1.0, range(m), range(m)), matrix(0.0, (m,1))
    ffdelays = graph.get_ffdelays()
    if type == 'Polynomial':
        coefs = graph.get_coefs()
        if not SO: coefs = coefs * spdiag([1.0/(j+2) for j in range(coefs.size[1])])
        parameters = matrix([[ffdelays], [coefs]])
        G = ue.objective_poly
    if type == 'Hyperbolic':
        ks = graph.get_ks()
        parameters = matrix([[ffdelays-div(ks[:,0],ks[:,1])], [ks]])
        G = ue.objective_hyper
    def F(x=None, z=None):
        if x is None: return 0, solver_init(U,r,random)
        if z is None:
            f, Df = G(P*x, z, parameters, 1)
            return f, Df*P
        f, Df, H = G(P*x, z, parameters, 1)
        return f, Df*P, P.T*H*P    
    x = solvers.cp(F, G=A, h=b, A=U, b=r)['x']
    if update:
        print 'Update link flows, delays in Graph.'; graph.update_linkflows_linkdelays(P*x)
        print 'Update path delays in Graph.'; graph.update_pathdelays()
        print 'Update path flows in Graph object.'; graph.update_pathflows(x)
    return x


def feasible_pathflows(graph, l_obs, obs=None, update=False, eq_constraints=None):
    """Attempts to find feasible pathflows given partial of full linkflows
    
    Parameters:
    ----------
    graph: Graph object
    l_obs: observations
    obs: indices of the observed links
    update: if True, update path flows in graph
    """
    P, n = linkpath_incidence(graph), graph.numpaths
    if eq_constraints is None: U, r = simplex(graph)
    else: U, r = eq_constraints
    if obs is not None: P2 = P[obs,:]
    C, d, q = spmatrix(-1.0, range(n), range(n)), matrix(0.0, (n,1)), -P2.trans()*l_obs
    x = solvers.qp(P2.trans()*P2, q, C, d, U, r)['x']
    if update:
        print 'Update link flows, delays in Graph.'; graph.update_linkflows_linkdelays(P*x)
        print 'Update path delays in Graph.'; graph.update_pathdelays()
        print 'Update path flows in Graph object.'; graph.update_pathflows(x)
    return x
    