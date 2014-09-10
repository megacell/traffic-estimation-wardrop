'''
Created on Sep 9, 2014

@author: jeromethai
'''

import numpy as np
import ue_solver as ue
import inverse_opt as invopt
from cvxopt import matrix, spdiag, solvers


def compute_delays(l, ffdelays, coefs):
    """Compute delays given linkflows l"""
    n, d = coefs.size
    delays = matrix(0.0, (n,1))
    for i in range(n):
        delays[i] = ffdelays[i] + coefs[i,:] * matrix(np.power(l[i],range(1,d+1)))
    return delays


def compute_cost(data, toll=None, SO=False):
    """Compute (tolled) UE equilibrium or SO equilibrium
    
    Return value
    ------------
    total cost
    linkflows
    """
    Aeq, beq, ffdelays, coefs, type = data
    if toll is None: toll = 0.0
    x = ue.solver(data=(Aeq, beq, ffdelays+toll, coefs, type), SO=SO)
    delays = compute_delays(x, ffdelays, coefs)
    return (delays.T*x)[0], x


def ty_solver(data, l, w_so, w_toll):
    """Solves the block (t,y) of the toll pricing model
    
    Parameters
    ----------
    data: Aeq, beq, ffdelays, coefs
    l: linkflows
    w_so: weight placed on the SO objective
    w_toll: weight placed on the min toll collection objective
    """
    Aeq, beq, ffdelays, coefs = data
    delays = compute_delays(l, ffdelays, coefs)
    n = len(l)
    p = Aeq.size[1]/n
    m = Aeq.size[0]/p
    c = matrix([(1.0+w_so+w_toll)*l, -beq])
    I = spdiag([1.0]*n)
    G1 = matrix([-I]*(p+1))
    G2 = matrix([Aeq.T, matrix(0.0, (n,p*m))])
    G = matrix([[G1],[G2]])
    h = [delays]*p
    h.append(matrix(0.0, (n,1)))
    h = matrix(h)
    return solvers.lp(c,G,h)['x'][range(n)]


def solver(data, w_so, w_toll, max_iter=5):
    """Solves the toll pricing problem
    
    Parameters
    ----------
    data: Aeq, beq, ffdelays, coefs
    w_so: weight on the SO objective
    w_toll: weight on the min toll collection objective
    max_iter: maximum number of iterations
    """
    Aeq, beq, ffdelays, coefs = data
    n = len(ffdelays)
    p = Aeq.size[1]/n
    toll = matrix(0.0, (n,1))
    for k in range(max_iter):
        l = invopt.x_solver(ffdelays+(1.0+w_toll/(w_so+1))*toll, coefs, Aeq, beq)
        print (compute_delays(l, ffdelays, coefs).T*l)[0]
        toll = ty_solver(data, l, w_so, w_toll)
    return toll


def main_solver(graph, theta, ws=[1e-8, 1e-6, 1e-4, 1e-2, 1.0], max_iter=5):
    """Main solver for the toll pricing model
    
    Parameters
    ----------
    graph: Graph object
    theta: coefficients of polynomial graph
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


if __name__ == '__main__':
    pass