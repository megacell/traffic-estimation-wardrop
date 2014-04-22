'''
Created on Apr 20, 2014

@author: jeromethai
'''

from cvxopt import matrix, spmatrix
import scipy.linalg as sla

def constraints(graph):
    """Construct constraints for the UE"""
    keys = {graph.links.keys()[i]: i for i in range(graph.numlinks)}
    entries = []; I = []; J = []; beq=[]
    for id1,node in graph.nodes.items():
        for id2,link in node.inlinks.items():
            entries.append(1); I.append(id1-1); J.append(keys[id2])
        for id2,link in node.outlinks.items():
            entries.append(-1); I.append(id1-1); J.append(keys[id2])
        beq.append( sum([od.flow for od in node.endODs.values()]) - sum([od.flow for od in node.startODs.values()]) )
    Aeq = spmatrix(entries,I,J)
    ind = findBasis(matrix(Aeq).trans())
    n = graph.numlinks
    A, b, Aeq, beq = spmatrix(-1.0, range(n), range(n)), matrix(0.0, (n,1)), Aeq[ind,:], matrix(beq)[ind]
    return keys, A, b, Aeq, beq


def findBasis(M):
    """Find the indices of the columns of M that form a basis or range(M)"""
    p,l,u = sla.lu(M)
    ind = [i for i in range(u.shape[0]) if u[i,i] != 0.0]
    if u[i,i] == 0:
        for j in range(i+1,u.shape[1]):
            if u[i,j] != 0.0: ind.append(j); break
    return ind


def solver(graph, update=False):
    """Solve for the UE and update graph if update==True"""
    
    keys, A, b, Aeq, beq = constraints(graph)
    
    type = graph.links.values()[0].delayfunc.type    
    
    if type == 'Affine':
        from cvxopt.solvers import qp
        entries = []; I=[]; q = []
        for id,link in graph.links.items():
            entries.append(link.delayfunc.slope); I.append(keys[id]); q.append(link.delayfunc.ffdelay)
        P = spmatrix(entries,I,I)
        sol = qp(P, matrix(q), A, b, Aeq, beq)['x']
    
    if update == True:
        for id,link in graph.links.items(): flow = sol[keys[id]]; link.flow, link.delay = flow, link.delayfunc.computeDelay(flow)
        for path in graph.paths.values():
            delay = 0.0
            for link in path.links: delay += link.delay
            path.delay = delay
        
    return keys, sol