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


def constraints(graph, linkflows_obs=None, indlinks_obs=None, soft=None):
    """Construct constraints for the UE link flow
    
    Parameters
    ----------
    graph: graph object
    linkflows_obs: vector of observed link flows (must be in the same order)
    indlinks_obs: list of indices of observed link flows (must be in the same order)
    soft: if provided, constraints to reconcile x^obs are switched to soft constraints
    
    Return value
    ------------
    Aeq: matrix of incidence nodes-links
    beq: matrix of OD flows at each node
    """
    m, n = graph.numnodes, graph.numlinks
    entries, I, J, beq = [], [], [], matrix(0.0, (m,1))
    
    for id1,node in graph.nodes.items():
        for id2,link in node.inlinks.items(): entries.append(1.0); I.append(id1-1); J.append(graph.indlinks[id2])
        for id2,link in node.outlinks.items(): entries.append(-1.0); I.append(id1-1); J.append(graph.indlinks[id2])
        beq[id1-1] = sum([od.flow for od in node.endODs.values()]) - sum([od.flow for od in node.startODs.values()])
    Aeq = spmatrix(entries, I, J, (m,n))
    
    num_obs = 0
    if linkflows_obs is not None and indlinks_obs is not None and soft is None:
        print 'Include observed link flows as hard constraints'
        num_obs = len(indlinks_obs)
        entries, I, J = np.ones(num_obs), range(num_obs), []
        for id in indlinks_obs: J.append(graph.indlinks[id])
        Aeq = sparse([Aeq, spmatrix(entries, I, J, (num_obs,n))])
        beq = matrix([beq, matrix(linkflows_obs)])
        
    M = matrix(Aeq); r = rank(M)
    if r < m+num_obs: print 'Remove {} redundant constraint(s)'.format(m+num_obs-r); ind = find_basis(M.trans())
    return Aeq[ind,:], beq[ind]


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
    return C, ind


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


def solver(graph, update=False, Aeq=None, beq=None, linkflows_obs=None, indlinks_obs=None, soft=None):
    """Find the UE link flow
    
    Parameters
    ----------
    graph: graph object
    update: if update==True: update link flows and link,path delays in graph
    Aeq: matrix of incidence nodes-links
    beq: matrix of OD flows at each node
    linkflows_obs: vector of observed link flows (must be in the same order)
    indlinks_obs: list of indices of observed link flows (must be in the same order)
    soft: if provided, constraints to reconcile x^obs are switched to soft constraints
    """
    
    if Aeq is None or beq is None: Aeq, beq = constraints(graph, linkflows_obs, indlinks_obs, soft)
    n = graph.numlinks
    A, b = spmatrix(-1.0, range(n), range(n)), matrix(0.0, (n,1))
    type = graph.links.values()[0].delayfunc.type    
        
    if type == 'Polynomial':
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
            if x is None: return 0, matrix(1.0, (n,1))
            #if min(x) <= 0.0: return None #implicit constraints
            f, Df = 0.0, matrix(0.0, (1,n))
            tmp2 = matrix(0.0, (n,1))
            for id,link in graph.links.items():
                i = graph.indlinks[id]
                tmp1 = matrix(np.power(x[i],range(degree+2)))
                f += ffdelays[i]*x[i] + coefs_i[i,:] * tmp1[2:degree+2]
                Df[i] = ffdelays[i] + coefs[i,:] * tmp1[1:degree+1]
                tmp2[i] = coefs_d[i,:] * tmp1[:degree]
            if linkflows_obs is not None and indlinks_obs is not None and soft is not None:
                obs = [graph.indlinks[id] for id in indlinks_obs]
                num_obs = len(obs)
                f += 0.5*soft*np.power(np.linalg.norm(x[obs]-linkflows_obs),2)
                I, J = [0]*num_obs, obs
                Df += soft*(spmatrix(x[obs],I,J, (1,n)) - spmatrix(linkflows_obs,I,J, (1,n)))
                tmp2 += spmatrix([soft]*num_obs,J,I, (n,1))
            if z is None: return f, Df
            H = spdiag(z[0] * tmp2)
            return f, Df, H
        
        if linkflows_obs is not None and indlinks_obs is not None and soft is not None:
            print 'Include observed link flows as soft constraints'
        linkflows = solvers.cp(F, G=A, h=b, A=Aeq, b=beq)['x']
    else:
        print 'Not implemented yet for non-polynomial delay functions'
        return
    
    if update:
        print 'Update link flows and link delays in Graph object.'; graph.update_linkflows_linkdlays(linkflows)
        print 'Update path delays in Graph object.'; graph.update_pathdelays()
        
    return linkflows
