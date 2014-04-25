'''
Created on Apr 20, 2014

@author: jeromethai
'''

from cvxopt import matrix, spmatrix, solvers
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
    M = matrix(Aeq); ind = range(m); r = rank(M)
    if r < m: print 'Remove {} redundant constraint(s)'.format(m-r); ind = find_basis(M.trans())
    return Aeq[ind,:], matrix(beq)[ind]


def solver(graph, update=True):
    """Find the UE link flow
    if update==True: update link flows and link,path delays in graph"""
    
    Aeq, beq = constraints(graph)
    n = graph.numlinks
    A, b = spmatrix(-1.0, range(n), range(n)), matrix(0.0, (n,1))
    
    type = graph.links.values()[0].delayfunc.type    
    
    if type == 'Affine':
        entries = []; I=[]; q = matrix(0.0, (graph.numlinks,1))
        for id,link in graph.links.items():
            entries.append(link.delayfunc.slope); I.append(graph.indlinks[id]); q[graph.indlinks[id]] = link.delayfunc.ffdelay
        P = spmatrix(entries,I,I)
        linkflows = solvers.qp(P, matrix(q), A, b, Aeq, beq)['x']
        
    if type == 'Other':
        pass
    
    if update == True:
        print 'Update link flows and link delays in Graph object.'
        for id,link in graph.links.items(): flow = linkflows[graph.indlinks[id]]; link.flow, link.delay = flow, link.delayfunc.compute_delay(flow)
        print 'Update path delays in Graph object.'
        for path in graph.paths.values():
            delay = 0.0
            for link in path.links: delay += link.delay
            path.delay = delay
        
    return linkflows


def unused_paths(graph, tol=1e-3):
    """Find unused paths given UE link delays in graph"""
    pathids = []
    for od in graph.ODs.values():
        mindelay = min([path.delay for path in od.paths.values()])
        [pathids.append((path.o, path.d, path.route)) for path in od.paths.values() if path.delay > mindelay*(1 + tol)]
    return pathids