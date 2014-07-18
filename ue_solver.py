'''
Created on Apr 20, 2014

@author: jeromethai
'''

import numpy as np
import scipy.sparse as sps
import scipy.io as sio
from cvxopt import matrix, spmatrix, solvers, spdiag, sparse
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
    
    degree = graph.links.values()[0].delayfunc.degree
    coefs, coefs_i, coefs_d = matrix(0.0, (n, degree)), matrix(0.0, (n, degree)), matrix(0.0, (n, degree))
    ffdelays = matrix(0.0, (n,1))
    for id,link in graph.links.items():
        i = graph.indlinks[id]
        ffdelays[i] = link.delayfunc.ffdelay
        for j in range(degree):
            coef = link.delayfunc.coef[j]
            coefs[i,j], coefs_i[i,j], coefs_d[i,j] = coef, coef/(j+2), coef*(j+1)
                    
    def F(x=None, z=None):
        if x is None: return 0, matrix(1.0/p, (p*n,1))
        y = matrix(0.0, (n,1))
        for k in range(p): y += x[k*n:(k+1)*n]
        f, Df = 0.0, matrix(0.0, (1,n))
        tmp2 = matrix(0.0, (n,1))
        for id,link in graph.links.items():
            i = graph.indlinks[id]
            tmp1 = matrix(np.power(y[i],range(degree+2)))
            f += ffdelays[i]*y[i] + coefs_i[i,:] * tmp1[2:degree+2]
            Df[i] = ffdelays[i] + coefs[i,:] * tmp1[1:degree+1]
            tmp2[i] = coefs_d[i,:] * tmp1[:degree]
        Df = matrix([[Df]]*p)
        if z is None: return f, Df
        H = spdiag(z[0] * tmp2)
        return f, Df, matrix([[H]*p]*p)
    
    x = solvers.cp(F, G=A, h=b, A=Aeq, b=beq)['x']
    linkflows = matrix(0.0, (n,1))
    for k in range(p): linkflows += x[k*n:(k+1)*n]
    
    if update:
        print 'Update link flows and link delays in Graph object.'; graph.update_linkflows_linkdelays(linkflows)
        print 'Update path delays in Graph object.'; graph.update_pathdelays()
    
    if full: return linkflows, x    
    return linkflows
