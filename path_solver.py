'''
Created on Jul 23, 2014

@author: jeromethai
'''

REG_EPS = 1e-9
TOL = 1e-5

import ue_solver as ue
from cvxopt import matrix, spmatrix, spdiag, solvers, div
import numpy.random as ra
import numpy as np
import logging
if logging.getLogger().getEffectiveLevel() >= logging.DEBUG:
    solvers.options['show_progress'] = False
else:
    solvers.options['show_progress'] = True
import Waypoints as WP
import rank_nullspace as rn
from util import find_basis
from kktsolver import get_kktsolver



def linkpath_incidence(graph):
    """Returns matrix of incidence link-path
    """
    I, J, entries = [], [], []
    for id1,link in graph.links.items():
        if len(link.paths) > 0:
            for id2 in link.paths.keys():
                I.append(graph.indlinks[id1])
                J.append(graph.indpaths[id2])
                entries.append(1.0)
        else:
            I.append(graph.indlinks[id1])
            J.append(0)
            entries.append(0.0)
    return spmatrix(entries, I, J)


def path_to_OD_simplex(graph):
    """Construct constraints for feasible path flows
    
    Return value
    ------------
    U: matrix of simplex constraints
    r: matrix of OD flows
    """
    I, J, r = [], [], matrix(0.0, (graph.numODs,1))
    for id1,od in graph.ODs.items():
        r[graph.indods[id1]] = od.flow
        for id2,path in od.paths.items():
            I.append(graph.indods[id1])
            J.append(graph.indpaths[id2])
    return spmatrix(1.0, I, J), r


def solver_init(U,r, random=False):
    """Initialize with a feasible point
    
    Parameters:
    ---------
    U: matrix of simplex constraints
    r: matrix of OD flows
    random: if true
    """
    n,m = U.size
    x0 = matrix(0.0, (m,1))
    for i in range(n):
        k = int(sum(U[i,:]))
        if random:
            tmp = ra.rand(k)
            tmp = r[i]*(tmp/sum(tmp))
        else: tmp = [r[i]/k]*k
        k = 0
        for j in range(m):
            if U[i,j] > 0.0: x0[j] = tmp[k]; k += 1
    return x0


def solver(graph, update=False, data=None, SO=False, random=False):
    """Solve for the UE equilibrium using link-path formulation
    
    Parameters
    ----------
    graph: graph object
    update: if True, update link and path flows in graph
    data: (P,U,r) 
            P: link-path incidence matrix
            U,r: simplex constraints 
    SO: if True compute SO
    random: if True, initialize with a random feasible point
    """
    type = graph.links.values()[0].delayfunc.type
    if data is None:
        P = linkpath_incidence(graph)
        U,r = path_to_OD_simplex(graph)
    else: P,U,r = data
    m = graph.numpaths
    A, b = spmatrix(-1.0, range(m), range(m)), matrix(0.0, (m,1))
    ffdelays = graph.get_ffdelays()
    if type == 'Polynomial':
        coefs = graph.get_coefs()
        if not SO: coefs = coefs * spdiag([1.0/(j+2) for j in range(coefs.size[1])])
        parameters = matrix([[ffdelays], [coefs]])
        G = ue.objective_poly
    if type == 'Hyperbolic':
        ks = graph.get_ks()
        parameters = matrix([[ffdelays-div(ks[:,0],ks[:,1])], [ks]])
        G = ue.objective_hyper
    def F(x=None, z=None):
        if x is None: return 0, solver_init(U,r,random)
        if z is None:
            f, Df = G(P*x, z, parameters, 1)
            return f, Df*P
        f, Df, H = G(P*x, z, parameters, 1)
        return f, Df*P, P.T*H*P    
    failed = True
    while failed:
        x = solvers.cp(F, G=A, h=b, A=U, b=r)['x']
        l = ue.solver(graph, SO=SO)
        error = np.linalg.norm(P * x - l,1) / np.linalg.norm(l,1)
        if error > TOL:
            print 'error={} > {}, re-compute path_flow'.format(error, TOL)
        else: failed = False
    if update:
        logging.info('Update link flows, delays in Graph.'); graph.update_linkflows_linkdelays(P*x)
        logging.info('Update path delays in Graph.'); graph.update_pathdelays()
        logging.info('Update path flows in Graph object.'); graph.update_pathflows(x)
    # assert if x is a valid UE/SO
    return x, l


def feasible_pathflows(graph, l_obs, obs=None, update=False,
                       with_cell_paths=False, with_ODs=False, x_true=None, wp_trajs=None):
    """Attempts to find feasible pathflows given partial of full linkflows
    
    Parameters:
    ----------
    graph: Graph object
    l_obs: observations of link flows
    obs: indices of the observed links
    update: if True, update path flows in graph
    with_cell_paths: if True, include cell paths as constraints
    with_ODs: if True, include ODs in the constraints if no with_cell_paths or in the objective if with_cell_paths
    """
    assert with_cell_paths or with_ODs # we must have some measurements!
    n = graph.numpaths
    # route to links flow constraints
    A, b = linkpath_incidence(graph), l_obs
    if obs: A = A[obs,:] # trim matrix if we have partial observations
    Aineq, bineq = spmatrix(-1.0, range(n), range(n)), matrix(0.0, (n,1)) # positive constraints
    if not with_cell_paths: # if just with ODs flow measurements:
        Aeq, beq = path_to_OD_simplex(graph) # route to OD flow constraints
    else: # if we have cellpath flow measurements:
        assert wp_trajs is not None
        Aeq, beq = WP.simplex(graph, wp_trajs) # route to cellpath flow constraints
        if with_ODs: # if we have ODs + cellpaths measurements
          T, d = path_to_OD_simplex(graph) # route to OD flow constraints included in objective
          A, b = matrix([A, T]), matrix([b, d]) # add the constraints to the objective
        
    if x_true is not None:
        err1 =  np.linalg.norm(A * x_true - b, 1) / np.linalg.norm(b, 1)
        err2 = np.linalg.norm(Aeq * x_true - beq) / np.linalg.norm(beq, 1)
        assert err1 < TOL, 'Ax!=b'
        assert err2 < TOL, 'Aeq x!=beq'
    # construct objective for cvxopt.solvers.qp
    Q, c = A.trans()*A, -A.trans()*b
    #x = solvers.qp(Q + REG_EPS*spmatrix(1.0, range(n), range(n)), c, Aineq, bineq, Aeq, beq)['x']
    # try with cvxopt.solvers.cp
    def qp_objective(x=None, z=None):
      if x is None: return 0, matrix(1.0, (n, 1))
      f = 0.5 * x.trans()*Q*x + c.trans() * x
      Df = (Q*x + c).trans()
      if z is None: return f, Df
      return f, Df, z[0]*Q
  
    dims = {'l': n, 'q': [], 's': []}
    x = solvers.cp(qp_objective, G=Aineq, h=bineq, A=Aeq, b=beq, 
        kktsolver=get_kktsolver(Aineq, dims, Aeq, qp_objective))['x']
    
    if update:
        logging.info('Update link flows, delays in Graph.'); graph.update_linkflows_linkdelays(P*x)
        logging.info('Update path delays in Graph.'); graph.update_pathdelays()
        logging.info('Update path flows in Graph object.'); graph.update_pathflows(x)
    return x, rn.rank(matrix([A, Aeq])), n
    