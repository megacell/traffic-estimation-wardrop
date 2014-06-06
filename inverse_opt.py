'''
Created on Jun 6, 2014

@author: jeromethai
'''

import numpy as np

def constraints(list_graphs, list_linkflows, degree):
    """Construct constraints for the Inverse Optimization
    (only available for polynomial delays)
    delay at link i is D_i(x) = ffdelays[i] + sum^{degree}_{k=1} theta[k-1]*(slope[i]*x)^k
    estimate the coefficients theta
    
    Parameters
    ----------
    list_graphs: list of graphs with same geometry, different OD demands
                 delay functions of links must contain ffdelays and slopes 
    
    list_linkflows: list of link flow vectors in equilibrium
    
    degree: degree of the polynomial function to estimate
    
    Return value
    ------------
    c, A, b, Aeq, beq: such that 
        min c'*x 
        s.t. A*x <= b
        Aeq*x = beq
    """
    
    if len(list_graphs) != len(list_linkflows): print 'ERROR: Lists not same length'; return
    
    graph = list_graphs[0]
    type = graph.links.values()[0].delayfunc.type 
    
    if type != 'Polynomial':
        print 'Inverse Optimization only available for polynomial delay functions'
        return
    
    N = len(list_graphs)
    m, n = graph.numnodes, graph.numlinks
    
    ffdelays, slopes = matrix(0.0, (n,1)), matrix(0.0, (n,1))
    for id,link in graph.links.items():
        i = graph.indlinks[id]
        ffdelays[i], slopes[i] = link.delayfunc.ffdelay, link.delayfunc.slope
    
    if type == 'Polynomial':
        tmp1, tmp2 = matrix(0.0, (degree,1)), matrix(0.0, (n,degree))
        tmp3, tmp4 = matrix(0.0, (n*N,m*N)), matrix(0.0, (m*N,1))
        for j, graph, linkflows in zip(range(N), list_graphs, list_linkflows):
            tmp1 += 1
        
        """
        for id,link in graph.links.items():
                i = graph.indlinks[id]
                tmp1 = matrix(np.power(x[i],range(degree+2)))
                f += ffdelays[i]*x[i] + coefs_int[i,:] * tmp1[range(2,degree+2)]
                Df[i] = ffdelays[i] + coefs[i,:] * tmp1[range(1,degree+1)]
                tmp2[i] = coefs_der[i,:] * tmp1[range(degree)]
        """
        
        
        return c, A, b, Aeq, beq
    


def solver(graph):
    """Solves the inverse optimization problem
    (only available for polynomial delays)
    """
    type = graph.links.values()[0].delayfunc.type
    
    if type != 'Polynomial':
        print 'Inverse Optimization only available for polynomial delay functions'
        return
    
    if type == 'Polynomial':
        return