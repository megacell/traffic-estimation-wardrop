'''
Created on Jun 6, 2014

@author: jeromethai
'''

import ue_solver as ue
import numpy as np
from cvxopt import matrix, spmatrix, solvers, spdiag, mul
from util import bisection


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


def solver(list_graphs, list_linkflows, degree, smooth, data=None, fvalue=False, full=False):
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
    fvalue: if True, returns value of the objective without regularization
    full: if True, returns the full x
    """
    type = list_graphs[0].links.values()[0].delayfunc.type
    N = len(list_graphs)
    
    if type != 'Polynomial':
        print 'Inverse Optimization only available for polynomial delay functions'
        return
    
    if type == 'Polynomial':
        if data is None: data = constraints(list_graphs, list_linkflows, degree)
        c, A, b = data
        P = spmatrix(smooth, range(degree), range(degree), (len(c),len(c)))
        sol = solvers.qp(P, c, G=A, h=b)
        theta = sol['x'][range(degree)]
        x = theta
        if full: x = sol['x']
        if fvalue:
            return x, sol['primal objective'] + (b.T*matrix(list_linkflows))[0] - 0.5*(theta.T*P[:degree,:degree]*theta)[0]
        
    return x


def x_solver(graph, y, Aeq, beq, linkflows_obs, indlinks_obs, soft=None):
    """Solves the single-level optimization problem in 'x'
    
    Parameters
    ----------
    graph: graph object
    y: dual variables
    Aeq: matrix of incidence nodes-links
    beq: matrix of OD flows at each node
    linkflows_obs: vector of observed link flows (must be in the same order)
    indlinks_obs: list of indices of observed link flows (must be in the same order)
    soft: reconcile x^obs are switched to soft constraints
    """
    n = graph.numlinks
    m = Aeq.size[0]
    if soft is None: m -= len(indlinks_obs)
    A, b = spmatrix(-1.0, range(n), range(n)), matrix(0.0, (n,1))
    type = graph.links.values()[0].delayfunc.type
    
    if type != 'Polynomial':
        print 'Not implemented yet'
        return
    
    if type == 'Polynomial':
        degree = graph.links.values()[0].delayfunc.degree
        ffdelays = matrix(0.0, (n,1))
        coefs, coefs1, coefs2 = matrix(0.0, (n, degree)), matrix(0.0, (n, degree)), matrix(0.0, (n, degree))
        for id,link in graph.links.items():
            i = graph.indlinks[id]
            ffdelays[i] = link.delayfunc.ffdelay
            for j in range(degree):
                coef = link.delayfunc.coef[j]
                coefs[i,j] = coef
                coefs1[i,j] = coef*(j+2)
                coefs2[i,j] = coef*(j+2)*(j+1)
        
        def F(x=None, z=None):
            if x is None: return 0, matrix(1.0, (n,1))
            f, Df = 0.0, matrix(0.0, (1,n))
            tmp2 = matrix(0.0, (n,1))
            for id,link in graph.links.items():
                i = graph.indlinks[id]
                tmp1 = matrix(np.power(x[i],range(degree+2)))
                f += ffdelays[i]*x[i] + coefs[i,:] * tmp1[range(2,degree+2)]
                Df[i] = ffdelays[i] + coefs1[i,:] * tmp1[range(1,degree+1)]
                tmp2[i] = coefs2[i,:] * tmp1[range(degree)]
            if soft is not None:
                obs = [graph.indlinks[id] for id in indlinks_obs]
                num_obs = len(obs)
                f += 0.5*soft*np.power(np.linalg.norm(x[obs]-linkflows_obs),2)
                I, J = [0]*num_obs, obs
                Df += soft*(spmatrix(x[obs],I,J, (1,n)) - spmatrix(linkflows_obs,I,J, (1,n)))
                tmp2 += spmatrix([soft]*num_obs,J,I, (n,1))
            if z is None: return f, Df
            H = spdiag(z[0] * tmp2)
            return f, Df, H
        
        #print Aeq.size, y.size
        #print m,n
        tmp3 = Aeq[:m,:n].T*y
        for id,link in graph.links.items():
            i = graph.indlinks[id]
            if tmp3[i] > ffdelays[i]:
                def G(x): return np.dot(coef[i,:], np.power(x, range(1,degree+1)))
                b[i] = -bisection(G, tmp3[i]-ffdelays[i], 0.0, 5.0*link.delayfunc.slope)
        
        linkflows = solvers.cp(F, G=A, h=b, A=Aeq, b=beq)['x']
        
    return linkflows


def solver_mis(list_graphs, list_linkflows_obs, indlinks_obs, degree, smooth, soft=None, max_iter=2, fvalue=False, alt=False):
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
    fvalue: if True, returns value of the objective without regularization
    """
    N, graph = len(list_graphs), list_graphs[0]
    n = graph.numlinks
    obs = [graph.indlinks[id] for id in indlinks_obs]
    #theta_init = matrix(np.ones(degree))/float(degree)
    theta_init = matrix(np.zeros(degree))
    theta_init[0] = 1.0
    beqs = []
    for j in range(N):
        tmp1, tmp2 = ue.constraints(list_graphs[j], list_linkflows_obs[j], indlinks_obs, soft)
        if j<1: Aeq = tmp1 # same node-link incidence matrix
        beqs.append(tmp2) # different demands
    
    m = Aeq.size[0]
    if soft is None: m -= len(indlinks_obs)
    
    slopes, ffdelays = matrix(0.0, (n,1)), matrix(0.0, (n,1))
    for id,link in graph.links.items():
        slopes[graph.indlinks[id]] = link.delayfunc.slope
        ffdelays[graph.indlinks[id]] = link.delayfunc.ffdelay # same slopes
    
    theta, f1, f2 = theta_init, [], []
    ys = [matrix(0.0, (m,1)) for j in range(N)]
    
    for k in range(max_iter):
        
        if k%2 == 0:
            for id, link in graph.links.items():
                i = graph.indlinks[id]
                link.delayfunc.coef = [ffdelays[i]*a*b for a,b in zip(theta, np.power(slopes[i], range(1,degree+1)))]
            if alt:
                list_linkflows = [x_solver(graph, ys[j], Aeq, beqs[j], list_linkflows_obs[j], indlinks_obs, soft) for j in range(N)]
            else:
                list_linkflows = [ue.solver(graph, False, Aeq, beqs[j], list_linkflows_obs[j], indlinks_obs, soft) for j in range(N)]
            if fvalue:
                error = 0
                for i in range(N):
                    e = list_linkflows[i][obs] - list_linkflows_obs[i]
                    error += e.T*e
                    f1.append(error[0])
        else:
            if fvalue:
                x, f = solver(list_graphs, list_linkflows, degree, smooth, fvalue=True, full=True)
                f2.append(f)
            else:
                x = solver(list_graphs, list_linkflows, degree, smooth, fvalue=False, full=True)
            theta = x[range(degree)]
            ys = [x[degree+j*m:degree+(j+1)*m] for j in range(N)]
            
    if fvalue: return theta, f1, f2
    return theta


def direct_solver(list_graphs, list_linkflows_obs, indlinks_obs, degree, smooth, soft=None):
    """Solves the inverse optimization problem with missing values
    
    Parameters
    ----------
    list_graphs: list of graphs with same geometry, different OD demands
                 delay functions of links must contain ffdelays and slopes 
    list_linkflows_obs: list of partially-observed link flow vectors in equilibrium
    indlinks_obs: indices of observed links
    degree: degree of the polynomial function to estimate
    smooth: regularization parameter on theta
    soft: regularization parameter for soft constraints
    """
    N, graph = len(list_graphs), list_graphs[0]
    n = graph.numlinks
    obs = [graph.indlinks[id] for id in indlinks_obs]
    theta_init = matrix(np.zeros(degree))
    theta_init[0] = 1.0
    beqs = []
    for j in range(N):
        tmp1, tmp2 = ue.constraints(list_graphs[j], list_linkflows_obs[j], indlinks_obs, soft)
        if j<1: Aeq = tmp1 # same node-link incidence matrix
        beqs.append(tmp2) # different demands
    
    m = Aeq.size[0]
    if soft is None: m -= len(indlinks_obs)
    C = Aeq[:m,:n]
    
    slopes, ffdelays = matrix(0.0, (n,1)), matrix(0.0, (n,1)) # get ff delays and slopes
    for id,link in graph.links.items():
        i = graph.indlinks[id]
        slopes[i] = link.delayfunc.slope
        ffdelays[i] = link.delayfunc.ffdelay # same slopes
        
    if type != 'Polynomial':
        print 'Not implemented yet'
        return
    
    if type == 'Polynomial':
        degree = graph.links.values()[0].delayfunc.degree
        x0 = matrix([matrix(1.0, (N*n, 1)), theta_init, matrix(0.0, (N*m, 1))])
        
        def F(x=None, z=None):
            if x is None: return N*n, x0
            delays, delays1 = matrix(0.0, (N*n, 1)), matrix(0.0, (N*n, 1))
            f, Df = matrix(0.0, (N*n+1, 1)), matrix(0.0, (N*n+1, N*n+degree+N*m))
            theta = x[N*n:N*n+degree]
            xs = [x[j*n:(j+1)*n] for j in range(N)]
            ys = [x[N*n+degree+j*m:N*n+degree+j*(m+1)] for j in range(N)]
            tmp2 = matrix(np.multiply(theta, range(1,degree+1)))
            tmp3 = matrix(np.multiply)
            
            for j in range(N):
                for id,link in graph.links.items():
                    i = graph.indlinks[id]
                    tmp1 = matrix(np.power(slopes[i]*xs[j][i], range(degree+1)))
                    delays[j*n+i] = ffdelays[i]*(1.0 + tmp1[1:].T*theta)
                    delays1[j*n+i] = slopes[i]*(tmp1[:-1].T*tmp2)
                    Df[0, j*n+i] += ffdelays[i]*(tmp1[1:].T*range(2,degree+2))
                Df[0, j*n:(j+1)*n] = ffdelays.T
                f[0] += delays[j*n:j*(n+1)].T*xs[j] - beqs[j].T*ys[j]
                f[j*n+1:(j+1)*n+1] = C.T*ys[j] - delays[j*n:(j+1)*n]
                Df[0, N*n+degree+j*m:N*n+degree+j*(m+1)] = -beqs[j]
                Df[j*n+1:(j+1)*n+1, j*n:(j+1)*n] = -spmatrix(delays1[j*n:(j+1)*n], range(n), range(n))
                Df[j*n+1:(j+1)*n+1, N*n+degree+j*m:N*n+degree+j*(m+1)] = C.T
                
            P = spmatrix(smooth, range(degree), range(degree))
            f[0] += 0.5*(theta.T*P*theta)
            Df[0, N*n:N*n+degree] = matrix(np.multiply(smooth,theta))
            
            if soft is not None:
                num_obs = len(obs)
                for j in range(N):
                    linkflows_obs = list_linkflows_obs[j]
                    f += 0.5*soft*np.power(np.linalg.norm(xs[j][obs]-linkflows_obs),2)
                    I, J = [0]*num_obs, obs
                    Df[0,:] += soft*(spmatrix(xs[j][obs],I,J, (1,n)) - spmatrix(linkflows_obs,I,J, (1,n)))
        
            if z is None: return f, Df
            return f, Df, spdiag(H*z)
        
        theta = solvers.cp(F, G=A, h=b, A=Aeq, b=beq)['x'][N*n:N*n+degree]
    
    return theta
