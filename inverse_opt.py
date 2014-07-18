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
    delay at link i is D_i(x) = ffdelays[i]*(1 + sum^d_k=1 theta[k-1]*(slope[i]*x)^k)
    estimate the coefficients theta
    
    Parameters
    ----------
    list_graphs: list of graphs with same geometry, different OD demands
                 delay functions of links must contain ffdelays and slopes 
    list_linkflows: list of link flow vectors in equilibrium
    degree: degree of the polynomial function to estimate
        
    Return value
    ------------
    c, A, b: such that min c'*x + r(x) s.t. A*x <= b
    """    
    graph = list_graphs[0]
    type, n = graph.links.values()[0].delayfunc.type, graph.numlinks 
    if type != 'Polynomial': print 'Delay functions must be polynomial'; return
    N = len(list_graphs)
    Aeq, beq = ue.constraints(graph)
    p = Aeq.size[1]/n
    m = Aeq.size[0]/p
    ffdelays, slopes = graph.get_ffdelays(), graph.get_slopes()     
    b = matrix([ffdelays]*p*N)
    tmp1, tmp2 = matrix(0.0, (degree, 1)), matrix(0.0, (p*N*n, degree))
    tmp3, tmp4 = matrix(0.0, (p*n*N, m*N*p)), matrix(0.0, (p*m*N, 1))
    for j, graph, linkflows in zip(range(N), list_graphs, list_linkflows): 
        tmp5 = mul(slopes, linkflows)
        for deg in range(degree):
            tmp6 = mul(tmp5**(deg+1), ffdelays)
            tmp1[deg] += (linkflows.T) * tmp6
            tmp2[j*n*p:(j+1)*n*p, deg] = -matrix([tmp6]*p)
        if j>0: Aeq, beq = ue.constraints(graph)
        tmp3[j*n*p:(j+1)*n*p, j*m*p:(j+1)*m*p] = Aeq.T
        tmp4[j*m*p:(j+1)*m*p, 0] = -beq
                    
    c, A = matrix([tmp1, tmp4]), matrix([[tmp2], [tmp3]])
    return c, A, b


def solver(list_graphs, list_linkflows, degree, smooth, full=False):
    """Solves the inverse optimization problem
    (only available for polynomial delays)
    
    Parameters
    ----------
    list_graphs: list of graphs with same geometry, different OD demands
                 delay functions of links must contain ffdelays and slopes 
    list_linkflows: list of link flow vectors in equilibrium
    degree: degree of the polynomial function to estimate
    smooth: regularization parameters on theta
    full: if True, returns the full x
    """
    type = list_graphs[0].links.values()[0].delayfunc.type
    if type != 'Polynomial': print 'Delay functions must be polynomial'; return
    c, A, b = constraints(list_graphs, list_linkflows, degree)
    P = spmatrix(smooth, range(degree), range(degree), (len(c),len(c)))
    sol = solvers.qp(P, c, G=A, h=b)
    if full: return sol['x']
    return sol['x'][range(degree)]


def x_solver(graph, Aeq, beq, soft, obs, l_obs):
    """
    optimization w.r.t. x_block
    
    Parameters
    ---------
    graph: graph object
    Aeq, beq: equality constraints of the ue program
    soft: parameter
    obs: indices of the observed links
    l_obs: observations
    """
    n = graph.numlinks
    p = Aeq.size[1]/n
    A, b = spmatrix(-1.0, range(p*n), range(p*n)), matrix(0.0, (p*n,1))
    ffdelays, coefs = graph.get_ffdelays(), graph.get_coefs()
    def F(x=None, z=None): return ue.objective(x, z, matrix([[ffdelays], [coefs]]), p)
    x = solvers.cp(F, G=A, h=b, A=Aeq, b=beq)['x']
    linkflows = matrix(0.0, (n,1))
    for k in range(p): linkflows += x[k*n:(k+1)*n]
    return linkflows


def x_solver_2(graph, y, Aeq, beq, linkflows_obs, indlinks_obs, soft):
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
            obs = [graph.indlinks[id] for id in indlinks_obs]
            num_obs = len(obs)
            f += 0.5*soft*np.power(np.linalg.norm(x[obs]-linkflows_obs),2)
            I, J = [0]*num_obs, obs
            Df += soft*(spmatrix(x[obs],I,J, (1,n)) - spmatrix(linkflows_obs,I,J, (1,n)))
            tmp2 += spmatrix([soft]*num_obs,J,I, (n,1))
            if z is None: return f, Df
            H = spdiag(z[0] * tmp2)
            return f, Df, H
        
        tmp3 = Aeq[:m,:n].T*y
        for id,link in graph.links.items():
            i = graph.indlinks[id]
            if tmp3[i] > ffdelays[i]:
                def G(x): return np.dot(coef[i,:], np.power(x, range(1,degree+1)))
                b[i] = -bisection(G, tmp3[i]-ffdelays[i], 0.0, 5.0*link.delayfunc.slope)
        
        linkflows = solvers.cp(F, G=A, h=b, A=Aeq, b=beq)['x']
        
    return linkflows


def x_solver_3(graph, y, Aeq, beq, linkflows_obs, indlinks_obs, soft):
    """Solves the single-level optimization problem in 'x'
    
    Parameters
    ----------
    graph: graph object
    y: dual variables
    Aeq, beq: constraints of the ue_solver
    linkflows_obs: vector of observed link flows (must be in the same order)
    indlinks_obs: list of indices of observed link flows (must be in the same order)
    soft: reconcile x^obs
    """
    n = graph.numlinks
    p = Aeq.size[1]/n
    m = Aeq.size[0]/p
    C = Aeq[:m, :n]
    A1, b1, b2 = spmatrix(-1.0, range(p*n), range(p*n)), matrix(0.0, (p*n,1)), matrix(0.0, (n,1))
    type = graph.links.values()[0].delayfunc.type
    
    ffdelays = matrix(0.0, (n,1))
    for id,link in graph.links.items(): ffdelays[graph.indlinks[id]] = link.delayfunc.ffdelay
    
    A2 = matrix([[spmatrix(-1.0, range(n), range(n))]]*p)
    #A = matrix([A1,A2])
    A = A2
    
    if type != 'Polynomial':
        print 'Not implemented yet'
        return
    
    if type == 'Polynomial':
        degree = graph.links.values()[0].delayfunc.degree
        coefs, coefs1, coefs2 = matrix(0.0, (n, degree)), matrix(0.0, (n, degree)), matrix(0.0, (n, degree))
        for id,link in graph.links.items():
            i = graph.indlinks[id]
            for j in range(degree):
                coef = link.delayfunc.coef[j]
                coefs[i,j], coefs1[i,j], coefs2[i,j] = coef, coef*(j+2), coef*(j+2)*(j+1)
        
        tmp = matrix(0.0, (n,p))
        for k in range(p): tmp[:,k] = C.T*y[k*m:(k+1)*m]
        for id,link in graph.links.items():
            i = graph.indlinks[id]
            tmp4 = max(tmp[i,:])
            if tmp4 > ffdelays[i]:
                def G(x): return np.dot(coef[i,:], np.power(x, range(1,degree+1)))
                b2[i] = -bisection(G, tmp4-ffdelays[i], 0.0, 5.0/link.delayfunc.slope)
        
        #b = matrix([b1,b2])
        b = b2
        
        def F(x=None, z=None):
            if x is None: return 0, matrix(1.0/p, (p*n,1))
            l = matrix(0.0, (n,1))
            for k in range(p): l += x[k*n:(k+1)*n]
            f, Df = 0.0, matrix(0.0, (1,n))
            tmp2 = matrix(0.0, (n,1))
            for id,link in graph.links.items():
                i = graph.indlinks[id]
                tmp1 = matrix(np.power(l[i],range(degree+2)))
                f += ffdelays[i]*l[i] + coefs[i,:] * tmp1[range(2,degree+2)]
                Df[i] = ffdelays[i] + coefs1[i,:] * tmp1[range(1,degree+1)]
                tmp2[i] = coefs2[i,:] * tmp1[range(degree)]
                
            obs = [graph.indlinks[id] for id in indlinks_obs]
            num_obs = len(obs)
            f += 0.5*soft*np.power(np.linalg.norm(l[obs]-linkflows_obs),2)
            I, J = [0]*num_obs, obs
            Df += soft*(spmatrix(l[obs],I,J, (1,n)) - spmatrix(linkflows_obs,I,J, (1,n)))
            tmp2 += spmatrix([soft]*num_obs,J,I, (n,1))
            Df = matrix([[Df]]*p)
            if z is None: return f, Df
            H = spdiag(z[0] * tmp2)
            return f, Df, matrix([[H]*p]*p)
        
        x = solvers.cp(F, G=A, h=b, A=Aeq, b=beq)['x']
        linkflows = matrix(0.0, (n,1))
        for k in range(p): linkflows += x[k*n:(k+1)*n]
        
    return linkflows





def solver_mis(list_graphs, list_linkflows_obs, indlinks_obs, degree, smooth, soft=1000.0, max_iter=10, fvalue=False):
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
    theta_init = matrix(np.zeros(degree))
    theta_init[0] = 1.0
    beqs = []
    for j in range(N):
        tmp1, tmp2 = ue.constraints(list_graphs[j])
        if j<1: Aeq = tmp1 # same node-link incidence matrix
        beqs.append(tmp2) # different demands
    
    p = Aeq.size[1]/n
    m = Aeq.size[0]/p
    
    slopes, ffdelays = matrix(0.0, (n,1)), matrix(0.0, (n,1))
    for id,link in graph.links.items():
        i = graph.indlinks[id]
        slopes[i] = link.delayfunc.slope
        ffdelays[i] = link.delayfunc.ffdelay # same slopes
    
    theta, f1, f2 = theta_init, [], []
    ys = [matrix(0.0, (m*p,1)) for j in range(N)]
    
    for k in range(max_iter):
        
        if k%2 == 0:
            for id, link in graph.links.items():
                i = graph.indlinks[id]
                link.delayfunc.coef = [ffdelays[i]*a*b for a,b in zip(theta, np.power(slopes[i], range(1,degree+1)))]
            #list_linkflows = [ue.solver(graph, ys[j], Aeq, beqs[j], list_linkflows_obs[j], indlinks_obs, soft) for j in range(N)]
            list_linkflows = [ue.solver(graph, False, Aeq, beqs[j], list_linkflows_obs[j], indlinks_obs, soft) for j in range(N)]

         
            if fvalue:
                for i in range(N):
                    e = list_linkflows[i][obs] - list_linkflows_obs[i]
                    f1.append((e.T*e)[0])
        else:
            if fvalue:
                x, f = solver(list_graphs, list_linkflows, degree, smooth, fvalue=True, full=True)
                f2.append(f)
            else:
                x = solver(list_graphs, list_linkflows, degree, smooth, fvalue=False, full=True)
            theta = x[range(degree)]
            ys = [x[degree+j*m*p:degree+(j+1)*m*p] for j in range(N)]
            
    if fvalue: return theta, f1, f2
    return theta
