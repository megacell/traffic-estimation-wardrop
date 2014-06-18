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
    c, A, b: such that 
        min c'*x + r(x)
        s.t. A*x <= b
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
        
        b = matrix(0.0, (n*N, 1))
        tmp1, tmp2 = matrix(0.0, (degree,1)), matrix(0.0, (N*n, degree))
        tmp3, tmp4 = matrix(0.0, (n*N, m*N)), matrix(0.0, (m*N,1))
        ind_start1, ind_end1 = [j*m for j in range(N)], [(j+1)*m for j in range(N)]
        ind_start2, ind_end2 = [j*n for j in range(N)], [(j+1)*n for j in range(N)]
        
        for j, graph, linkflows in zip(range(N), list_graphs, list_linkflows): 
            
            tmp5 = mul(slopes, linkflows)
            for deg in range(degree):
                tmp6 = mul(tmp5**(deg+1), ffdelays)
                tmp1[deg] += (linkflows.T) * tmp6
                tmp2[ind_start2[j]:ind_end2[j], deg] = spmatrix(-ffdelays,range(n),range(n))*tmp6
                
            if j>0: C, beq = ue.constraints(graph)
            tmp3[ind_start2[j]:ind_end2[j], ind_start1[j]:ind_end1[j]] = C.T
            tmp4[ind_start1[j]:ind_end1[j], 0] = -beq
            
            b[ind_start2[j]:ind_end2[j]] = ffdelays
            
        scale = (sum(abs(tmp4)) * float(len(tmp1))) / (sum(abs(tmp1)) * float(len(tmp4)))
        c, A = matrix([scale*tmp1, tmp4]), matrix([[scale*tmp2], [tmp3]])
    
        #for deg in range(degree): A[deg, deg] = -1.0
        return c, A, scale*b


def solver(list_graphs, list_linkflows, degree, smooth, data=None):
    """Solves the inverse optimization problem
    (only available for polynomial delays)
    
    Parameters
    ----------
    list_graphs: list of graphs with same geometry, different OD demands
                 delay functions of links must contain ffdelays and slopes 
    
    list_linkflows: list of link flow vectors in equilibrium
    
    degree: degree of the polynomial function to estimate
        
    smooth: regularization parameters on theta
        
    data: constraints and objective for the optimization problem
    """
    type = list_graphs[0].links.values()[0].delayfunc.type
    
    if type != 'Polynomial':
        print 'Inverse Optimization only available for polynomial delay functions'
        return
    
    if type == 'Polynomial':
        if data is None: data = constraints(list_graphs, list_linkflows, degree)
        c, A, b = data
        print c.size
        print A.size
        print b.size
        P = spmatrix(smooth, range(degree), range(degree), (len(c),len(c)))
        #x = solvers.lp(c, G=A, h=b)['x']
        x = solvers.qp(P, c, G=A, h=b)['x']
        
    return x[range(degree)]


def solver_mis(list_graphs, list_linkflows_obs, indlinks_obs, degree, smooth, soft=None, max_iter=2):
    """Solves the inverse optimization problem with missing values
    (only available for polynomial delays)
    
    Parameters
    ----------
    list_graphs: list of graphs with same geometry, different OD demands
                 delay functions of links must contain ffdelays and slopes 
    
    list_linkflows_obs: list of partially-observed link flow vectors in equilibrium
    
    indlinks_obs: indices of observed links
    
    degree: degree of the polynomial function to estimate
    
    smooth: regularization parameter on theta
    
    soft: regularization parameter for soft constraints
    
    max_iter: maximum number of iterations
    """
    N, graph = len(list_graphs), list_graphs[0]
    n = graph.numlinks
    #theta_init = matrix(np.ones(degree))/float(degree)
    theta_init = matrix(np.zeros(degree))
    theta_init[0] = 1.0
    beqs = []
    for j in range(N):
        tmp1, tmp2 = ue.constraints(list_graphs[j], list_linkflows_obs[j], indlinks_obs, soft)
        if j<1: Aeq = tmp1 # same node-link incidence matrix
        beqs.append(tmp2) # different demands
    
    slopes, ffdelays = matrix(0.0, (n,1)), matrix(0.0, (n,1))
    for id,link in graph.links.items():
        slopes[graph.indlinks[id]] = link.delayfunc.slope
        ffdelays[graph.indlinks[id]] = link.delayfunc.ffdelay # same slopes
    
    theta = theta_init
    
    for k in range(max_iter):
        
        if k%2 == 0:
            for id, link in graph.links.items():
                i = graph.indlinks[id]
                link.delayfunc.coef = [ffdelays[i]*a*b for a,b in zip(theta, np.power(slopes[i], range(1,degree+1)))]
            list_linkflows = [ue.solver(graph, False, Aeq, beqs[j], list_linkflows_obs[j], indlinks_obs, soft) for j in range(N)]
        else:
            theta = solver(list_graphs, list_linkflows, degree, smooth)
        
    return theta