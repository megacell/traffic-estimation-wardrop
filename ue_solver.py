'''
Created on Apr 20, 2014

@author: jeromethai
'''

from cvxopt import matrix, spmatrix
import scipy.linalg as sla
from rank_nullspace import rank


def index_links(graph):
    """Index links of the graph with integers from 0 to numlinks-1"""
    return {graph.links.keys()[i]: i for i in range(graph.numlinks)}


def constraints(graph):
    """Construct constraints for the UE link flow
    
     Return value
    ------------
    indlinks: dict of link id to index
    A: matrix -eye(numlinks)
    b: matrix zeros(numlinks,1)
    Aeq: matrix of incidence nodes-links
    beq: matrix of incidence nodes-ODs
    """
    indlinks = index_links(graph)
    entries, I, J, beq = [], [], [], []
    for id1,node in graph.nodes.items():
        for id2,link in node.inlinks.items(): entries.append(1.0); I.append(id1-1); J.append(indlinks[id2])
        for id2,link in node.outlinks.items(): entries.append(-1.0); I.append(id1-1); J.append(indlinks[id2])
        beq.append( sum([od.flow for od in node.endODs.values()]) - sum([od.flow for od in node.startODs.values()]) )
    Aeq = spmatrix(entries,I,J)
    M = matrix(Aeq); m = graph.numnodes; ind = range(m); r = rank(M)
    if r < m: print 'Remove {} redundant constraint(s)'.format(m-r); ind = find_basis(M.trans())
    n = graph.numlinks
    A, b, Aeq, beq = spmatrix(-1.0, range(n), range(n)), matrix(0.0, (n,1)), Aeq[ind,:], matrix(beq)[ind]
    return indlinks, A, b, Aeq, beq


def find_basis(M):
    """Find the indices of the columns of M that form a basis or range(M)"""
    p,l,u = sla.lu(M)
    ind = [i for i in range(u.shape[0]) if u[i,i] != 0.0]
    if u[i,i] == 0:
        for j in range(i+1,u.shape[1]):
            if u[i,j] != 0.0: ind.append(j); break
    return ind


def solver(graph, update=False):
    """Find the UE link flow
    if update==True: update link flows and link,path delays in graph"""
    
    indlinks, A, b, Aeq, beq = constraints(graph)
    
    type = graph.links.values()[0].delayfunc.type    
    
    if type == 'Affine':
        from cvxopt.solvers import qp
        entries = []; I=[]; q = []
        for id,link in graph.links.items():
            entries.append(link.delayfunc.slope); I.append(indlinks[id]); q.append(link.delayfunc.ffdelay)
        P = spmatrix(entries,I,I)
        linkflows = qp(P, matrix(q), A, b, Aeq, beq)['x']
        
    if type == 'Other':
        pass
    
    if update == True:
        for id,link in graph.links.items(): flow = linkflows[indlinks[id]]; link.flow, link.delay = flow, link.delayfunc.compute_delay(flow)
        for path in graph.paths.values():
            delay = 0.0
            for link in path.links: delay += link.delay
            path.delay = delay
        
    return indlinks, linkflows