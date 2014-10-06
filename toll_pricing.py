'''
Created on Sep 9, 2014

@author: jeromethai
'''

import numpy as np
import ue_solver as ue
import inverse_opt as invopt
from cvxopt import matrix, spdiag, solvers, div


def compute_delays(l, ffdelays, pm, type='Polynomial'):
    """Compute delays given linkflows l"""
    n, d = pm.size
    delays = matrix(0.0, (n,1))
    if type == 'Polynomial':
        for i in range(n):
            delays[i] = ffdelays[i] + pm[i,:] * matrix(np.power(l[i],range(1,d+1)))
    if type == 'Hyperbolic':
        ks = matrix([[ffdelays-div(pm[:,0],pm[:,1])], [pm]])
        for i in range(n):
            delays[i] = ks[i,0] + ks[i,1]/(ks[i,2]-l[i])
    return delays


def compute_cost(data, toll=None, SO=False):
    """Compute:
    - untolled UE
    - untolled SO
    - tolled UE
    
    Return value
    ------------
    total cost
    linkflows
    """
    Aeq, beq, ffdelays, pm, type = data
    if toll is None: toll = 0.0
    x = ue.solver(data=(Aeq, beq, ffdelays+toll, pm, type), SO=SO)
    delays = compute_delays(x, ffdelays, pm, type)
    return (delays.T*x)[0], x


def ty_solver(data, l, w_toll, full=False):
    """Solves the block (t,y) of the toll pricing model
    
    Parameters
    ----------
    data: Aeq, beq, ffdelays, coefs
    l: linkflows
    w_toll: weight on the toll collected
    full: if True, return the whole solution
    """
    Aeq, beq, ffdelays, coefs = data
    delays = compute_delays(l, ffdelays, coefs)
    n = len(l)
    p = Aeq.size[1]/n
    m = Aeq.size[0]/p
    c = matrix([(1.0+w_toll)*l, -beq])
    I = spdiag([1.0]*n)
    G1 = matrix([-I]*(p+1))
    G2 = matrix([Aeq.T, matrix(0.0, (n,p*m))])
    G = matrix([[G1],[G2]])
    h = [delays]*p
    h.append(matrix(0.0, (n,1)))
    h = matrix(h)
    if full: return solvers.lp(c,G,h)['x']
    return solvers.lp(c,G,h)['x'][range(n)]


def solver(data, w_so, w_toll, max_iter=5, full=False):
    """Solves the toll pricing problem
    
    Parameters
    ----------
    data: Aeq, beq, ffdelays, coefs
    w_so: weight on the SO objective
    w_toll: weight on the toll collected
    max_iter: maximum number of iterations
    full: if False, return toll, if True, return toll, y, l
    """
    Aeq, beq, ffdelays, coefs = data
    n = len(ffdelays)
    p = Aeq.size[1]/n
    toll = matrix(100.0, (n,1))
    for k in range(max_iter):
        l = invopt.x_solver(ffdelays+((1.+w_toll)/(1.+w_so))*toll, coefs, Aeq, beq)
        print (compute_delays(l, ffdelays, coefs).T*l)[0]
        toll = ty_solver(data, l, w_toll)
    if full:
        sol = ty_solver(data, l, w_toll, True)
        return sol[:n], sol[n:], l # return toll, y, l
    return toll


def main_solver(graph, theta, ws=[1e-6], max_iter=5):
    """Main solver for the toll pricing model
    
    Parameters
    ----------
    graph: Graph object
    theta: coefficients of polynomial link delays
    ws: list of weights on the toll pricing (suppose w_toll=w_so)
    max_iter: maximum number of iterations
    """
    ffdelays, slopes = graph.get_ffdelays(), graph.get_slopes()
    Aeq, beq = ue.constraints(graph)
    coefs = invopt.compute_coefs(ffdelays, slopes, theta)
    data = (Aeq, beq, ffdelays, coefs)
    k, min_cost, toll_collected = len(ws), np.inf, 0.0
    for i in range(k):
        t = solver(data, ws[i], ws[i], max_iter)
        cost, x = compute_cost((Aeq, beq, ffdelays, coefs, 'Polynomial'), t)
        if cost < min_cost:
            toll = t
            min_cost = cost
            toll_collected = (toll.T*x)[0]
            weight = ws[i]
    return toll, min_cost, toll_collected, weight


def compute_gap(l, y, ks, toll, beq, scale):
    """Compute the gap function
    
    Parameters
    ----------
    l: link flow
    y: estimated dual variables
    ks: matrix of size (n,degree) where ks[i,:] are the coefs of the polynomial delay on link i
    beq: ue constraint Aeq*x=beq 
    toll: vector of tolls
    scale: scaling factor for the gap function
    """
    n, d = ks.size
    print toll.size
    print l.size
    print beq.size
    print y.size
    gap = toll.T*l - beq.T*y
    for i in range(n):
        tmp = matrix(np.power(l[i],range(1,d+1)))
        gap += ks[i,:] * tmp
    return scale*gap[0]


def multi_objective_solver(graph, theta, ws, max_iter=5):
    """Multi-objective solver for the toll pricing model
    
    Parameters
    ----------
    graph: graph object
    theta: coefficients of polynomial link delays
    ws: list of weights for w_so and w_toll
    max_iter: maximum number of iterations
    """
    ffdelays, slopes = graph.get_ffdelays(), graph.get_slopes()
    Aeq, beq = ue.constraints(graph)
    coefs = invopt.compute_coefs(ffdelays, slopes, theta)
    ue_cost = compute_cost((Aeq, beq, ffdelays, coefs, 'Polynomial'))[0]
    so_cost = compute_cost((Aeq, beq, ffdelays, coefs, 'Polynomial'), 0.0, True)[0]
    data = (Aeq, beq, ffdelays, coefs)
    r_gap = matrix(0.0, (len(ws), len(ws)))
    loss_est = matrix(0.0, (len(ws), len(ws)))
    toll_est = matrix(0.0, (len(ws), len(ws)))
    loss_res = matrix(0.0, (len(ws), len(ws)))
    toll_res = matrix(0.0, (len(ws), len(ws)))
    for i,w_so in enumerate(ws):
        for j,w_toll in enumerate(ws):
            t,y,l = solver(data, w_so, w_toll, max_iter, True)
            r_gap[i,j] = compute_gap(l, y, matrix([[ffdelays], [coefs]]), t, beq, 1/ue_cost)
            toll_est[i,j] = (t.T*l)[0]
            delays = compute_delays(l, ffdelays, coefs, 'Polynomial')
            cost = (delays.T*l)[0]
            loss_est[i,j] = (cost-so_cost)/(ue_cost-so_cost)
            cost, x = compute_cost((Aeq, beq, ffdelays, coefs, 'Polynomial'), t)
            toll_res[i,j] = (t.T*x)[0]
            loss_res[i,j] = (cost-so_cost)/(ue_cost-so_cost)
    return r_gap, toll_est, loss_est, toll_res, loss_res


if __name__ == '__main__':
    pass