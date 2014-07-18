'''
Created on Apr 20, 2014

@author: jeromethai
'''

import numpy as np
from cvxopt import matrix, spmatrix, solvers, spdiag, mul
from rank_nullspace import rank
from util import find_basis


def constraints(graph):
    """Construct constraints for the UE link flow
    
    Parameters
    ----------
    graph: graph object
    
    Return value
    ------------
    Aeq, beq: equality constraints Aeq*x = beq
    """
    C, ind = get_nodelink_incidence(graph)
    ds = [get_demands(graph, ind, id) for id,node in graph.nodes.items() if len(node.endODs) > 0]
    p = len(ds)
    m,n = C.size
    Aeq, beq = spmatrix([], [], [], (p*m,p*n)), matrix(ds)
    for k in range(p): Aeq[k*m:(k+1)*m, k*n:(k+1)*n] = C
    return Aeq, beq


def get_nodelink_incidence(graph):
    """
    get node-link incidence matrix
    
    Parameters
    ----------
    graph: graph object
    
    Return value
    ------------
    C: matrix of incidence node-link
    ind: indices of a basis formed by the rows of C
    """
    m, n = graph.numnodes, graph.numlinks
    entries, I, J = [], [], []
    for id1,node in graph.nodes.items():
        for id2,link in node.inlinks.items(): entries.append(1.0); I.append(id1-1); J.append(graph.indlinks[id2])
        for id2,link in node.outlinks.items(): entries.append(-1.0); I.append(id1-1); J.append(graph.indlinks[id2])
    C = spmatrix(entries, I, J, (m,n))
    M = matrix(C); r = rank(M)
    if r < m: print 'Remove {} redundant constraint(s)'.format(m-r); ind = find_basis(M.trans())
    return C[ind,:], ind


def get_demands(graph, ind, node_id):
    """
    get demands for all OD pairs sharing the same destination
    
    Parameters
    ----------
    graph: graph object
    ind: indices of a basis formed by the rows of C
    node_id: id of the destination node
    """
    d = matrix(0.0, (graph.numnodes,1))
    for OD in graph.nodes[node_id].endODs.values():
        d[node_id-1] += OD.flow
        d[OD.o-1] = -OD.flow
    return d[ind]


def objective(x, z, coefs, p, soft=0.0, obs=None, l_obs=None):
    """Objective function of UE program
    f(x) = sum_i f_i(l_i) (+ 0.5*soft*||l[obs]-l_obs||^2)
    f_i(u) = sum_{k=1}^degree coefs[i,k] u^k
    with l = sum_w x_w
    
    Parameters
    ----------
    x,z: variables for the F(x,z) function for cvxopt.solvers.cp
    coefs: matrix of size (n,degree) 
    p: number of w's
    soft: parameter
    obs: indices of the observed links
    l_obs: observations
    """
    n, d = coefs.size
    if x is None: return 0, matrix(1.0/p, (p*n,1))
    l = matrix(0.0, (n,1))
    for k in range(p): l += x[k*n:(k+1)*n]
    f, Df, H = 0.0, matrix(0.0, (1,n)), matrix(0.0, (n,1))
    for i in range(n):
        tmp = matrix(np.power(l[i],range(d+1)))
        f += coefs[i,:] * tmp[1:]
        Df[i] = coefs[i,:] * mul(tmp[:-1], matrix(range(1,d+1)))
        H[i] = coefs[i,1:] * mul(tmp[:-2], matrix(range(2,d+1)), matrix(range(1,d)))
    Df = matrix([[Df]]*p)
    
    if soft != 0.0:
        num_obs, e = len(obs), l[obs]-l_obs
        f += 0.5*soft*e.T*e
        Df += soft*spmatrix(e, [0]*num_obs, obs, (1,n))
        H += spmatrix([soft]*num_obs, obs, [0]*num_obs, (n,1))
    
    if z is None: return f, Df
    return f, Df, matrix([[spdiag(z[0] * H)]*p]*p)


def solver(graph, update=False, Aeq=None, beq=None, full=False):
    """Find the UE link flow
    
    Parameters
    ----------
    graph: graph object
    update: if update==True: update link flows and link,path delays in graph
    Aeq: matrix of incidence nodes-links
    beq: matrix of OD flows at each node
    full: if full=True, also return x (link flows per OD pair)
    """
    if Aeq is None or beq is None: Aeq, beq = constraints(graph)
    n = graph.numlinks
    p = Aeq.size[1]/n
    A, b = spmatrix(-1.0, range(p*n), range(p*n)), matrix(0.0, (p*n,1))
    type = graph.links.values()[0].delayfunc.type
    if type != 'Polynomial': print 'Delay functions must be polynomial'; return
    ffdelays, coefs = graph.get_ffdelays(), graph.get_coefs()
    coefs_i = coefs * spdiag([1.0/(j+2) for j in range(coefs.size[1])])
    def F(x=None, z=None): return objective(x, z, matrix([[ffdelays], [coefs_i]]), p)
    x = solvers.cp(F, G=A, h=b, A=Aeq, b=beq)['x']
    linkflows = matrix(0.0, (n,1))
    for k in range(p): linkflows += x[k*n:(k+1)*n]
    
    if update:
        print 'Update link flows, delays in Graph.'; graph.update_linkflows_linkdelays(linkflows)
        print 'Update path delays in Graph.'; graph.update_pathdelays()
    
    if full: return linkflows, x    
    return linkflows
