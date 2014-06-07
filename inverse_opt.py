'''
Created on Jun 6, 2014

@author: jeromethai
'''

import ue_solver as ue
import numpy as np
from cvxopt import matrix, spmatrix, solvers, spdiag, mul


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
    c, A, b, C, beq: such that 
        min c'*x 
        s.t. A*x <= b
        C*x = beq
    """
    
    if len(list_graphs) != len(list_linkflows): print 'ERROR: Lists not same length'; return
    
    graph = list_graphs[0]
    type = graph.links.values()[0].delayfunc.type 
    
    if type != 'Polynomial':
        print 'Inverse Optimization only available for polynomial delay functions'
        return
    
    N = len(list_graphs)
    C, beq = ue.constraints(graph)
    m, n = C.size
    
    ffdelays, slopes = matrix(0.0, (n,1)), matrix(0.0, (n,1))
    for id,link in graph.links.items():
        i = graph.indlinks[id]
        ffdelays[i], slopes[i] = link.delayfunc.ffdelay, link.delayfunc.slope
    
    if type == 'Polynomial':
        
        b, Aeq = matrix(0.0, (n*N + degree, 1)), matrix(0.0, (1, N*m + degree))
        tmp1, tmp2 = matrix(0.0, (degree,1)), matrix(0.0, (N*n + degree, degree))
        tmp3, tmp4 = matrix(0.0, (n*N + degree, m*N)), matrix(0.0, (m*N,1))
        ind_start1, ind_end1 = [j*m for j in range(N)], [(j+1)*m for j in range(N)]
        ind_start2, ind_end2 = [j*n for j in range(N)], [(j+1)*n for j in range(N)]
        
        for j, graph, linkflows in zip(range(N), list_graphs, list_linkflows): 
            
            tmp5 = mul(slopes, linkflows)
            for deg in range(degree):
                tmp6 = tmp5**(deg+1)
                tmp1[deg] += (linkflows.T) * tmp6
                tmp2[degree+ind_start2[j]:degree+ind_end2[j], deg] = -tmp6
                
            if j>0: C, beq = ue.constraints(graph)
            tmp3[degree+ind_start2[j]:degree+ind_end2[j], ind_start1[j]:ind_end1[j]] = C.T
            tmp4[ind_start1[j]:ind_end1[j], 0] = -beq
            
            b[degree+ind_start2[j]:degree+ind_end2[j]] = ffdelays
            
        scale = (sum(abs(tmp4)) * float(len(tmp1))) / (sum(abs(tmp1)) * float(len(tmp4)))
        c, A, beq = matrix([scale*tmp1, tmp4]), matrix([[scale*tmp2], [tmp3]]), matrix([1.0])
    
        for deg in range(degree):
            Aeq[deg] = 1.0
            A[deg, deg] = -1.0
    
        return c, A, scale*b, Aeq, beq


def solver(list_graphs, list_linkflows, degree, data=None):
    """Solves the inverse optimization problem
    (only available for polynomial delays)
    """
    type = list_graphs[0].links.values()[0].delayfunc.type
    
    if type != 'Polynomial':
        print 'Inverse Optimization only available for polynomial delay functions'
        return
    
    if type == 'Polynomial':
        if data is None: data = constraints(list_graphs, list_linkflows, degree)
        c, A, b, Aeq, beq = data
        x = solvers.lp(c, G=A, h=b, A=Aeq, b=beq)['x']
        
    return x[range(degree)]