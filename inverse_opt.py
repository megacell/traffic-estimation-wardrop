'''
Created on Jun 6, 2014

@author: jeromethai
'''
import ue_solver as ue
import numpy as np
from cvxopt import matrix, spmatrix, solvers, mul, spdiag
from multiprocessing import Pool
import shortest_paths as sh
import rank_nullspace as rn


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
        Aeq, tmp = ue.constraints(graph, True)
        beqs.append(tmp)
    return Aeq, beqs, ffdelays, slopes


def constraints(data, ls, degree):
    """Construct constraints for the Inverse Optimization with polynomial delay
    delay at link i is D_i(x) = ffdelays[i]*(1 + sum^d_k=1 theta[k-1]*(slope[i]*x)^k)
    estimate the coefficients theta
    
    Parameters
    ----------
    data: data for the inverse optimization problem, given by get_data
    ls: list of link flow vectors in equilibrium
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
    for j, beq, linkflows in zip(range(N), beqs, ls): 
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


def solver(data, ls, degree, full=False):
    """Solves the inverse optimization problem
    
    Parameters
    ----------
    data: Aeq, beqs, ffdelays, slopes from get_data(graphs)
    ls: list of link flow vectors in equilibrium
    degree: degree of the polynomial function to estimate
    full: if False, just return theta, if True, return the whole primal solution
    """
    c, A, b = constraints(data, ls, degree)
    print A.size
    r = rn.rank(matrix(A))
    print r
    if full: return solvers.lp(c, G=A, h=b)['x']
    return solvers.lp(c, G=A, h=b)['x'][range(degree)]


def x_solver(ffdelays, coefs, Aeq, beq, w_obs=0.0, obs=None, l_obs=None, w_gap=1.0):
    """
    optimization w.r.t. x_block:
    min F(x)'x + r(x) s.t. x in K
    
    Parameters
    ----------
    ffdelays: matrix of freeflow delays from graph.get_ffdelays()
    coefs: matrix of coefficients from graph.get_coefs()
    Aeq, beq: equality constraints of the ue program
    w_obs: weight on the observation residual
    obs: indices of the observed links
    l_obs: observations
    w_gap: weight on the gap function
    """
    n = len(ffdelays)
    p = Aeq.size[1]/n    
    A, b = spmatrix(-1.0, range(p*n), range(p*n)), matrix(0.0, (p*n,1))
    def F(x=None, z=None): return ue.objective_poly(x, z, matrix([[ffdelays], [coefs]]), p, w_obs, obs, l_obs, w_gap)
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


def solver_mis(data, ls_obs, obs, degree, w_obs=1000.0, max_iter=3, full=False, w_gap=1.0):
    """Solves the inverse optimization problem with missing values
    
    Parameters
    ----------
    data: Aeq, beqs, ffdelays, slopes from get_data(graphs)
    ls_obs: list of partially-observed link flow vectors in equilibrium
    obs: indlinks ids of observed links
    degree: degree of the polynomial function to estimate
    w_obs: weight on the observation residual
    max_iter: maximum number of iterations
    full: if False, just return thetam if True, return theta, ys, ls
    w_gap: weight on the gap function
    """
    Aeq, beqs, ffdelays, slopes = data
    N, n = len(beqs), len(ffdelays)
    p = Aeq.size[1]/n
    m = Aeq.size[0]/p
    theta = matrix(np.zeros(degree)); theta[0] = 1.0 # initial theta
    for k in range(max_iter):
        coefs = compute_coefs(ffdelays, slopes, theta)
        ls = [x_solver(ffdelays, coefs, Aeq, beqs[j], w_obs, obs, ls_obs[j], w_gap) for j in range(N)]
        theta = solver(data, ls, degree)
    if full:
        sol = solver(data, ls, degree, True)
        return sol[range(degree)], [sol[j*m*p:((j+1)*m*p)] for j in range(N)], ls #return theta, ys, ls
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


def smm_helper(data):
    """Helper function for solver_mis_multi"""
    ffdelays, coefs, Aeq, beq, soft, obs, l_obs = data
    return x_solver(ffdelays, coefs, Aeq, beq, soft, obs, l_obs)


def solver_mis_multi_thread(data, ls_obs, obs, degree, soft=1000.0, max_iter=3, processes=1):
    """Solves the inverse optimization problem with missing values using several processes
    
    Parameters
    ----------
    data: Aeq, beqs, ffdelays, slopes from get_data(graphs)
    ls_obs: list of partially-observed link flow vectors in equilibrium
    obs: indlinks ids of observed links
    degree: degree of the polynomial function to estimate
    soft: regularization parameter for soft constraints
    max_iter: maximum number of iterations
    processes: number of processes
    """
    Aeq, beqs, ffdelays, slopes = data
    N, n = len(beqs), len(ffdelays)
    p = Aeq.size[1]/n
    m = Aeq.size[0]/p
    theta = matrix(np.zeros(degree)); theta[0] = 1.0 #initial theta
    for k in range(max_iter):
        coefs = compute_coefs(ffdelays, slopes, theta)
        pool = Pool(processes=processes)
        ls = pool.map(smm_helper, zip([ffdelays]*N, [coefs]*N, [Aeq]*N, beqs, [soft]*N, [obs]*N, ls_obs))
        theta = solver(data, ls, degree)
    return theta


def main_solver(graphs, ls_obs, obs, degree, soft=1000.0, max_iter=3):
    """Solves solve_mis with the best parameter
    
    Parameters
    ----------
    graphs: list of graphs
    ls_obs: list of partially-observed link flow vectors in equilibrium
    obs: indlinks ids of observed links
    degree: degree of the polynomial function to estimate
    soft: regularization parameter for soft constraints
    max_iter: maximum number of iterations
    """
    data, n_obs = get_data(graphs), len(obs)
    Aeq, beqs, ffdelays, slopes = data
    type, N, min_e, n = 'Polynomial', len(beqs), np.inf, len(ffdelays)
    if n_obs < n: theta = solver_mis(data, ls_obs, obs, degree, soft, max_iter)
    #if n_obs < n: theta = solver_mis_multi_thread(data, ls_obs, obs, degree, i*np.ones(degree), soft, max_iter, 4)
    if n_obs == n: theta = solver(data, ls_obs, degree)
    coefs = compute_coefs(ffdelays, slopes, theta)
    xs = [ue.solver(data=(Aeq, beqs[j], ffdelays, coefs, type)) for j in range(N)]
    return theta, xs


def compute_scaling(graphs, ls_obs, obs, data):
    """Scale the objectives for multi-objective optimization
    
    Parameters
    ----------
    graphs: list of graphs
    ls_obs: list of observed links
    obs: indlinks ids of observed links
    data: obtained from get_data
    
    Return value
    ------------
    scale: vector of scaling factor such that
    scale[0]: scaling for the observation residual
    scale[1]: scaling for the gap function
    """
    N = len(graphs)
    scale = matrix(0.0, (2,1))
    Aeq, beqs, ffdelays, slopes = data
    theta = matrix([1.0]) #initial theta
    coefs = compute_coefs(ffdelays, slopes, theta)
    type = 'Polynomial'
    xs = [ue.solver(data=(Aeq, beqs[k], ffdelays, coefs, type)) for k in range(N)]
    scale[0] = 2.*np.linalg.norm(matrix(ls_obs)-matrix([x[obs] for x in xs]),2)**-2
    
    graph = graphs[0]
    total_delay = 0.0
    sinks = [id for id,node in graph.nodes.items() if len(node.endODs) > 0]
    for t in sinks:
        sources = [od[0] for od in graph.ODs.keys() if od[1]==t]
        As = sh.mainKSP(graph, sources, t, 1)
        for s,path in As.items():
            delay = 0.0
            for i in range(len(path[0])-1): delay += graph.links[(path[0][i],path[0][i+1],1)].delay
            demand = 0.0
            for g in graphs: demand += g.ODs[(s,t)].flow
            total_delay += 5.*delay*demand
    scale[1] = 1/total_delay
    return scale


def compute_gap(ls, ys, ks, beqs, scale):
    """Compute the gap function
    
    Parameters
    ----------
    ls: list of estimated link flows
    ys: list of estimated dual variables
    ks: matrix of size (n,degree) where ks[i,:] are the coefs of the polynomial delay on link i
    beqs: list of ue constraints Aeq*x=beq
    scale: scaling factor for the gap function
    """
    n, d = ks.size
    gap = 0.0
    for l,y,beq in zip(ls,ys,beqs):
        gap -= beq.T*y
        for i in range(n):
            tmp = matrix(np.power(l[i],range(1,d+1)))
            gap += ks[i,:] * tmp
    return scale*gap[0]
    


def multi_objective_solver(graphs, ls_obs, obs, degree, w_multi, max_iter=3):
    """Multi-objective optimization for the latency inference problem with the following weights:
    w1: weight on the observation residual
    w2: weight on the gap function
    w1 + w2 = 1
    
    Parameters
    ----------
    graphs: list of graphs
    ls_obs: list of partially-observed link flow vectors in equilibrium
    obs: indlinks ids of observed links
    degree: degree of the polynomial function to estimate
    w_multi: list of weights on the observation residual
    max_iter: maximum number of iterations
    """
    data = get_data(graphs)
    scale = compute_scaling(graphs, ls_obs, obs, data)
    w_obs = [scale[0]*w for w in w_multi] # weights on the observation residual
    w_gap = [scale[1]*(1-w) for w in w_multi] # weights on the gap function
    #softs = [(scale[0]*w)/(scale[1]*(1-w)) for w in w_multi]
    Aeq, beqs, ffdelays, slopes = data
    type, N = 'Polynomial', len(beqs)
    r_gap, r_obs, x_est, thetas = [], [], [], []
    for w1,w2 in zip(w_obs,w_gap):
        theta, ys, ls = solver_mis(data, ls_obs, obs, degree, 0.0*np.ones(degree), w1, max_iter, True, w2)
        coefs = compute_coefs(ffdelays, slopes, theta)
        r_gap.append(compute_gap(ls, ys, matrix([[ffdelays], [coefs]]), beqs, scale[1]))
        r_obs.append(scale[0]*0.5*np.linalg.norm(matrix(ls_obs)-matrix([l[obs] for l in ls]), 2)**2)
        xs = [ue.solver(data=(Aeq, beqs[k], ffdelays, coefs, type)) for k in range(N)]
        x_est.append(matrix(xs))
        thetas.append(theta)
    return r_gap, r_obs, x_est, thetas
