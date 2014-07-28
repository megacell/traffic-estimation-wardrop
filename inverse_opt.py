'''
Created on Jun 6, 2014

@author: jeromethai
'''
import ue_solver as ue
import numpy as np
from cvxopt import matrix, spmatrix, solvers, mul, spdiag


def get_data(graphs):
    """Get data for the inverse optimization problem
    
    Return value:
    Aeq: UE equality constraints for graphs[end]
    beqs: list of UE equality constraints for graph in graphs
    ffdelays: matrix of freeflow delays from graph.get_ffdelays()
    slopes: matrix of slopes from graph.get_slopes()
    """
    graph = graphs[0]
    ffdelays, slopes = graph.get_ffdelays(), graph.get_slopes()
    beqs = []
    for graph in graphs:
        Aeq, tmp = ue.constraints(graph)
        beqs.append(tmp)
    return Aeq, beqs, ffdelays, slopes


def constraints(data, flow_vectors, degree):
    """Construct constraints for the Inverse Optimization with polynomial delay
    delay at link i is D_i(x) = ffdelays[i]*(1 + sum^d_k=1 theta[k-1]*(slope[i]*x)^k)
    estimate the coefficients theta
    
    Parameters
    ----------
    data: data for the inverse optimization problem, given by get_data
    flow_vectors: list of link flow vectors in equilibrium
    degree: degree of the polynomial function to estimate
        
    Return value
    ------------
    c, A, b: such that min c'*x + r(x) s.t. A*x <= b
    """    
    Aeq, beqs, ffdelays, slopes = data
    N, n = len(beqs), len(ffdelays)
    p = Aeq.size[1]/n
    m = Aeq.size[0]/p
    b = matrix([ffdelays]*p*N)
    tmp1, tmp2 = matrix(0.0, (degree, 1)), matrix(0.0, (p*N*n, degree))
    tmp3, tmp4 = matrix(0.0, (p*n*N, m*N*p)), matrix(0.0, (p*m*N, 1))
    for j, beq, linkflows in zip(range(N), beqs, flow_vectors): 
        tmp5 = mul(slopes, linkflows)
        for deg in range(degree):
            tmp6 = mul(tmp5**(deg+1), ffdelays)
            tmp1[deg] += (linkflows.T) * tmp6
            tmp2[j*n*p:(j+1)*n*p, deg] = -matrix([tmp6]*p)
        tmp3[j*n*p:(j+1)*n*p, j*m*p:(j+1)*m*p] = Aeq.T
        tmp4[j*m*p:(j+1)*m*p, 0] = -beq
                    
    c, A = matrix([tmp1, tmp4]), matrix([[tmp2], [tmp3]])
    A = matrix([A, spmatrix(-1.0, range(degree), range(degree), (degree, degree+m*N*p))])
    b = matrix([b, matrix([0.0]*degree)])
    return c, A, b


def solver(graphs, flow_vectors, degree, smooth, data=None, full=False):
    """Solves the inverse optimization problem
    
    Parameters
    ----------
    graphs: list of graphs with same geometry, different OD demands
                 delay functions of links must contain ffdelays and slopes 
    flow_vectors: list of link flow vectors in equilibrium
    degree: degree of the polynomial function to estimate
    smooth: regularization parameters on theta
    full: if True, returns the full x
    """
    if data is None: data = get_data(graphs)
    c, A, b = constraints(data, flow_vectors, degree)
    P = spmatrix(smooth, range(degree), range(degree), (len(c),len(c)))
    sol = solvers.qp(P, c, G=A, h=b)
    if full: return sol['x']
    return sol['x'][range(degree)]


def x_solver(ffdelays, coefs, Aeq, beq, soft, obs, l_obs):
    """
    optimization w.r.t. x_block
    
    Parameters
    ----------
    ffdelays: matrix of freeflow delays from graph.get_ffdelays()
    coefs: matrix of coefficients from graph.get_coefs()
    Aeq, beq: equality constraints of the ue program
    soft: parameter
    obs: indices of the observed links
    l_obs: observations
    """
    n = len(ffdelays)
    p = Aeq.size[1]/n    
    A, b = spmatrix(-1.0, range(p*n), range(p*n)), matrix(0.0, (p*n,1))
    def F(x=None, z=None): return ue.objective_poly(x, z, matrix([[ffdelays], [coefs]]), p, soft, obs, l_obs)
    x = solvers.cp(F, G=A, h=b, A=Aeq, b=beq)['x']
    linkflows = matrix(0.0, (n,1))
    for k in range(p): linkflows += x[k*n:(k+1)*n]
    return linkflows


def xy_solver(ffdelays, coefs, Aeq, beq, soft, obs, l_obs, a):
    """
    optimization w.r.t. (x,y) block
    
    Parameters
    ----------
    ffdelays: matrix of freeflow delays from graph.get_ffdelays()
    coefs: matrix of coefficients from graph.get_coefs()
    Aeq, beq: equality constraints of the ue program
    soft: parameter
    obs: indices of the observed links
    l_obs: observations
    a: (nX1)-matrix s.t. D_i(u) ~ ffdelay[i] + a[i]*u for link i
    """
    n = len(ffdelays)
    p = Aeq.size[1]/n
    m = Aeq.size[0]/p
    A, b = spmatrix(-1.0, range(p*n), range(p*n)), matrix(0.0, (p*n,1))
    A = matrix([[A], [matrix(0.0, (p*n,p*m))]])
    parameters = matrix([[ffdelays], [coefs]])
    def F(x=None, z=None):
        if x is None: return 0, matrix([matrix(1.0/p, (p*n,1)), matrix(0.0, (p*m,1))])
        if z is None:
            f, Df = ue.objective_poly(x[:p*n], z, parameters, p, soft, obs, l_obs)
            return f + beq.T*x[p*n:], matrix([[Df], [beq.T]])
        f, Df, H = ue.objective_poly(x[:p*n], z, parameters, p, soft, obs, l_obs)
        return f + beq.T*x[p*n:], matrix([[Df], [beq.T]]), spdiag([H, matrix(0.0, (p*m,p*m))])
    A0, b0 = matrix([[matrix([[spdiag(-a)]*p]*p)], [Aeq.T]]), matrix([ffdelays]*p)
    A, b = matrix([A, A0]), matrix([b, b0])
    Aeq = matrix([[Aeq], [matrix(0.0, (p*m,p*m))]])
    x = solvers.cp(F, G=A, h=b, A=Aeq, b=beq)['x']
    linkflows = matrix(0.0, (n,1))
    for k in range(p): linkflows += x[k*n:(k+1)*n]
    return linkflows


def compute_coefs(ffdelays, slopes, theta):
    """Compute the matrix of coefficients
    
    Parameters
    ----------
    ffdelays: matrix of freeflow delays from graph.get_ffdelays()
    slopes: matrix of slopes from graph.get_slopes()
    theta: parameters
    """
    n, d = len(ffdelays), len(theta)
    coefs = matrix(0.0, (n,d))
    for i in range(n):
        coef = [ffdelays[i]*a*b for a,b in zip(theta, np.power(slopes[i], range(1,d+1)))]
        coefs[i,:] = matrix(coef).T
    return coefs


def solver_mis(graphs, ls_obs, obs, degree, smooth, soft=1000.0, max_iter=3):
    """Solves the inverse optimization problem with missing values
    
    Parameters
    ----------
    graphs: list of graphs with same geometry, different OD demands
    ls_obs: list of partially-observed link flow vectors in equilibrium
    obs: indlinks ids of observed links
    degree: degree of the polynomial function to estimate
    smooth: regularization parameter on theta
    soft: regularization parameter for soft constraints
    max_iter: maximum number of iterations
    """
    data = get_data(graphs)
    Aeq, beqs, ffdelays, slopes = data
    N, n = len(graphs), len(ffdelays)
    p = Aeq.size[1]/n
    m = Aeq.size[0]/p
    theta = matrix(np.zeros(degree)); theta[0] = 1.0
    for k in range(max_iter):
        coefs = compute_coefs(ffdelays, slopes, theta)
        ls = [x_solver(ffdelays, coefs, Aeq, beqs[j], soft, obs, ls_obs[j]) for j in range(N)]
        theta = solver(graphs, ls, degree, smooth, data)
    return theta

"""
alternative 1:
a = mul((theta.T*matrix(np.power(0.85, range(degree))))*slopes, ffdelays) #can replace by 1.0
ls = [xy_solver(ffdelays, coefs, Aeq, beqs[j], soft, obs, ls_obs[j], a) for j in range(N)]

alternative 2:
ls_old = ls; ls = []
for j in range(N):
    a = matrix(0.0, (n,1))
    for i in range(n):
        a[i] = coefs[i,:]*matrix(np.power(0.27*ls_old[j][i], range(degree)))
    ls.append(xy_solver(ffdelays, coefs, Aeq, beqs[j], soft, obs, ls_obs[j], a))
"""