'''
Created on Apr 20, 2014

@author: jeromethai
'''

import numpy as np
import scipy.sparse as sps
import scipy.io as sio
from cvxopt import matrix, spmatrix, solvers, spdiag
from rank_nullspace import rank
from util import find_basis


def constraints(graph):
    """Construct constraints for the UE link flow
    
    Return value
    ------------
    Aeq: matrix of incidence nodes-links
    beq: matrix of OD flows at each node
    """
    m = graph.numnodes;
    entries, I, J, beq = [], [], [], matrix(0.0, (m,1))
    
    for id1,node in graph.nodes.items():
        for id2,link in node.inlinks.items(): entries.append(1.0); I.append(id1-1); J.append(graph.indlinks[id2])
        for id2,link in node.outlinks.items(): entries.append(-1.0); I.append(id1-1); J.append(graph.indlinks[id2])
        beq[id1-1] = sum([od.flow for od in node.endODs.values()]) - sum([od.flow for od in node.startODs.values()])
    Aeq = spmatrix(entries,I,J)
    M = matrix(Aeq); r = rank(M)
    if r < m: print 'Remove {} redundant constraint(s)'.format(m-r); ind = find_basis(M.trans())
    return Aeq[ind,:], beq[ind]


def solver(graph, update=True, Aeq=None, beq=None, linkflows_obs=None, indlinks_obs=None):
    """Find the UE link flow
    
    Parameters
    ----------
    update: if update==True: update link flows and link,path delays in graph
    Aeq: matrix of incidence nodes-links
    beq: matrix of OD flows at each node
    linkflows_obs: vector of observed link flows (sorted)
    indlinks_obs: list of indices of observed link flows (sorted)
    """
    
    if Aeq is None or beq is None: Aeq, beq = constraints(graph)
    n = graph.numlinks
    A, b = spmatrix(-1.0, range(n), range(n)), matrix(0.0, (n,1))
    
    if linkflows_obs is not None and indlinks_obs is not None:
        print 'Include observed link flows in UE computation'
        num_obs = len(indlinks_obs)
        entries, I, J = np.ones(num_obs), range(num_obs), []
        for id in indlink_obs: J.append(graph.indlinks[id])
        Aeq = sparse([Aeq, spmatrix(entries,I,J)])
        beq = matrix([beq, linkflows_obs])
    
    type = graph.links.values()[0].delayfunc.type    
    
    if type == 'Affine':
        entries = []; I=[]; q = matrix(0.0, (graph.numlinks,1))
        for id,link in graph.links.items():
            entries.append(link.delayfunc.slope); I.append(graph.indlinks[id]); q[graph.indlinks[id]] = link.delayfunc.ffdelay
        P = spmatrix(entries,I,I)
        linkflows = solvers.qp(P, matrix(q), A, b, Aeq, beq)['x']
        
    if type == 'Polynomial':
        degree = graph.links.values()[0].delayfunc.degree
        coefs, coefs_int, coefs_der = matrix(0.0, (n, degree)), matrix(0.0, (n, degree)), matrix(0.0, (n, degree))
        ffdelays = matrix(0.0, (n,1))
        for id,link in graph.links.items():
            i = graph.indlinks[id]
            ffdelays[i] = link.delayfunc.ffdelay
            for j in range(degree):
                coef = link.delayfunc.coef[j]
                coefs[i,j] = coef
                coefs_int[i,j] = coef/(j+2)
                coefs_der[i,j] = coef*(j+1)
                        
        def F(x=None, z=None):
            if x is None: return 0, matrix(1.0, (n,1))
            #if min(x) <= 0.0: return None #implicit constraints
            f = 0.0
            Df = matrix(0.0, (1,n))
            tmp2 = matrix(0.0, (n,1))
            for id,link in graph.links.items():
                i = graph.indlinks[id]
                tmp1 = matrix(np.power(x[i],range(degree+2)))
                f += ffdelays[i]*x[i] + coefs_int[i,:] * tmp1[range(2,degree+2)]
                Df[i] = ffdelays[i] + coefs[i,:] * tmp1[range(1,degree+1)]
                tmp2[i] = coefs_der[i,:] * tmp1[range(degree)]
            if z is None: return f, Df
            H = spdiag(z[0] * tmp2)
            return f, Df, H
            
        linkflows = solvers.cp(F, G=A, h=b, A=Aeq, b=beq)['x']
        
                
    if type == 'Other':
        pass
    
    if update:
        print 'Update link flows and link delays in Graph object.'; graph.update_linkflows_linkdlays(linkflows)
        print 'Update path delays in Graph object.'; graph.update_pathdelays()
        
    return linkflows


def unused_paths(graph, tol=1e-3):
    """Find unused paths given UE link delays in graph"""
    pathids = []
    for od in graph.ODs.values():
        mindelay = min([path.delay for path in od.paths.values()])
        [pathids.append((path.o, path.d, path.route)) for path in od.paths.values() if path.delay > mindelay*(1 + tol)]
    return pathids


def save_mat(filepath, name, graph):
    """Save sparse matrices for the UE problem in solve_ue_data.mat file
    solve_ue_data.mat contains:
    C: incidence matrix node-link
    d: in-out-OD flows in each node
    If type_delayfunc = Affine:
    p: slopes of the delay function for each node
    q: constant coefficient
    """
    type = graph.links.values()[0].delayfunc.type
    entries, I, J, beq = [], [], [], np.zeros((graph.numnodes, 1)) 
    for id1,node in graph.nodes.items():
            for id2,link in node.inlinks.items(): entries.append(1.0); I.append(id1-1); J.append(graph.indlinks[id2])
            for id2,link in node.outlinks.items(): entries.append(-1.0); I.append(id1-1); J.append(graph.indlinks[id2])
            beq[id1-1] = sum([od.flow for od in node.endODs.values()]) - sum([od.flow for od in node.startODs.values()])
    Aeq = sps.coo_matrix((entries,(I,J)))
            
    if type == 'Affine':
        p, q = np.zeros((graph.numlinks, 1)), np.zeros((graph.numlinks, 1))
        for id,link in graph.links.items(): p[graph.indlinks[id]], q[graph.indlinks[id]] = link.delayfunc.slope, link.delayfunc.ffdelay
        sio.savemat(filepath + name + '.mat', mdict={'C': Aeq, 'd': beq, 'p': p, 'q': q})
        
    if type == 'Other':
        pass