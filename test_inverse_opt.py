'''
Created on Jun 6, 2014

@author: jeromethai
'''

import numpy as np
import ue_solver as ue
import inverse_opt as invopt
from test_graph import los_angeles
import matplotlib.pyplot as plt
from cvxopt import matrix
from util import add_noise

a, b = 3.5, 3.0
coef = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
degree = len(coef)
graph = los_angeles(coef, 'Polynomial')[0]

    
def display_results(true_linkflows, est_linkflows, true_theta, best_theta, delaytype):
    """Display results
    
    Parameters
    ----------
    true_linkflows: list of true linkflows
    est_linkflows: list of estimated linkflows
    true_theta: true parameters
    best_theta: best parameters
    """
    error = np.linalg.norm(matrix(true_linkflows) - matrix(est_linkflows))
    xdata = np.linspace(0.0, 2.5, num=100)
    vals = [1+(best_theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    if delaytype == 'Polynomial':
        true_vals = [1+(true_theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    if delaytype == 'Hyperbolic':
        a,b = true_theta
        true_vals = [1 - a/b + a/(b-x) for x in xdata]
    plt.plot(xdata, vals, 'r', label='estimate')
    plt.plot( xdata, true_vals, 'b', label='true')
    plt.xlabel('Link flow')
    plt.ylabel('Delay')
    plt.title(r'Estimated delay function. l2-norm error: {:.3f}'.format(error))
    plt.legend()
    plt.show()
    
    
def get_graphs_linkflows(theta, noisy=False, delaytype='Polynomial'):
    """Given parameters theta, get L.A. graphs and associated UE linkflows for polynomial delay functions
    """
    g1, g2, g3, g4 = los_angeles(theta, delaytype, noisy)
    l1, l2, l3, l4 = ue.solver(g1), ue.solver(g2), ue.solver(g3), ue.solver(g4)
    return g1, g2, g3, g4, l1, l2, l3, l4
    
    
def test1(delaytype):
    """find smooth parameter that minimizes ||x_true-x|| without noise and with all observations
    """
    if delaytype == 'Polynomial': true_theta = coef
    if delaytype == 'Hyperbolic': true_theta = (a,b)
    g1, g2, g3, g4, l1, l2, l3, l4 = get_graphs_linkflows(true_theta, delaytype=delaytype)
    min_error = np.inf
    for i in [0.01, 1.0, 10.0, 100.0]:
        theta = invopt.solver([g1, g2, g3, g4], [l1, l2, l3, l4], degree, i*np.ones(degree))
        g1, g2, g3, g4, x1, x2, x3, x4 = get_graphs_linkflows(theta)
        error = np.linalg.norm(matrix([l1,l2,l3,l4])-matrix([x1,x2,x3,x4]))
        if error < min_error:
            best_smooth, min_error, best_theta, y1, y2, y3, y4 = i, error, theta, x1, x2, x3, x4
    print best_smooth
    print best_theta
    display_results([l1, l2, l3, l4], [y1, y2, y3, y4], true_theta, best_theta, delaytype)
    

def test2(delaytype):
    """find smooth parameter using cross validation without noise and with all observations
    avg(smooth) = [27.5]*6 for poly delays, avg(smooth) = [2.7]*6 for hyper delay functions
    """
    if delaytype == 'Polynomial': true_theta = coef
    if delaytype == 'Hyperbolic': true_theta = (a,b)
    g1, g2, g3, g4, l1, l2, l3, l4 = get_graphs_linkflows(true_theta, delaytype=delaytype)
    gs, ls = [g1,g2,g3,g4], [l1,l2,l3,l4]
    best_smooth = [0.0]*4
    
    for index in range(4):
        min_error = np.inf
        for i in [0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]:
            theta = invopt.solver(gs[:index]+gs[index+1:], ls[:index]+ls[index+1:], degree, i*np.ones(degree))
            l = get_graphs_linkflows(theta)[4+index]
            error = np.linalg.norm(ls[index]-l)
            if error < min_error: best_smooth[index], min_error = i, error    
    print best_smooth


def test3(indlinks_obs):
    """Test x_solver, compute_lower, and compute_coefs
    """
    g1, g2, g3, g4, l1, l2, l3, l4 = get_graphs_linkflows(coef)
    obs = [g4.indlinks[id] for id in indlinks_obs]
    ffdelays, slopes, coefs = g4.get_ffdelays(), g4.get_slopes(), g4.get_coefs()
    Aeq, beq = ue.constraints(g4)
    C = ue.nodelink_incidence(g4)[0]
    lower = invopt.compute_lower(C, matrix(10.0, (Aeq.size[0],1)), ffdelays, slopes, coefs)
    print lower
    data = invopt.get_data([g1,g2,g3,g4])
    c, A, b = invopt.constraints(data, [l1, l2, l3, l4], degree)
    print c
    print invopt.x_solver(ffdelays, coefs, Aeq, beq, 1000.0, obs, l4[obs], lower)
    print invopt.compute_coefs(ffdelays, slopes, coef)
    
    
def test4(indlinks_obs, noisy, delaytype):
    """find parameters that minimizes the distance between x^obs_true in NOISY case
    and x^obs generated by each candidate function with PARTIAL observation
    """
    if delaytype == 'Polynomial': true_theta = coef
    if delaytype == 'Hyperbolic': true_theta = (a,b)
    g1, g2, g3, g4, l1, l2, l3, l4 = get_graphs_linkflows(true_theta, delaytype=delaytype)
    obs = [g1.indlinks[id] for id in indlinks_obs]
    z1, z2, z3, z4 = l1, l2, l3, l4
    if noisy:
        l1, l2, l3, l4 = add_noise(l1, 1/30.0), add_noise(l2, 1/30.0), add_noise(l3, 1/30.0), add_noise(l4, 1/30.0)
        g1, g2, g3, g4 = los_angeles(true_theta, 'Polynomial', True)
    min_error = np.inf
    for i in [0.01]:
        theta = invopt.solver_mis([g1, g2, g3, g4], [l1[obs], l2[obs], l3[obs], l4[obs]], 
                          indlinks_obs, degree, i*np.ones(degree))
        g1, g2, g3, g4, x1, x2, x3, x4 = get_graphs_linkflows(theta, noisy)
        e = np.linalg.norm(matrix([l1[obs],l2[obs],l3[obs],l4[obs]])-matrix([x1[obs],x2[obs],x3[obs],x4[obs]]))
        if e < min_error:
            best_smooth, min_error, best_theta = i, e, theta
            y1, y2, y3, y4 = x1, x2, x3, x4         
    print best_smooth
    print best_theta
    l1, l2, l3, l4 = z1, z2, z3, z4
    display_results([l1, l2, l3, l4], [y1, y2, y3, y4], true_theta, best_theta, delaytype)
    
    
def test5(indlinks_obs, faulty, noisy=False):
    """Cross validation to find a sensor that has been attacked
    results when obs = [(17,24),(24,40),(14,21),(16,23)] and (24,40) has been attacked
    """
    g1, g2, g3, g4, l1, l2, l3, l4 = get_graphs_linkflows(theta_true)
    if noisy:
        g1, g2, g3, g4 = los_angeles(theta_true, 'Polynomial', True)
        l1, l2, l3, l4 = add_noise(l1, 1/30.0), add_noise(l2, 1/30.0), add_noise(l3, 1/30.0), add_noise(l4, 1/30.0)
    faulty_id = g1.indlinks[faulty]
    l1[faulty_id], l2[faulty_id], l3[faulty_id], l4[faulty_id] = 0.5*l1[faulty_id], 0.5*l2[faulty_id], 0.5*l3[faulty_id], 0.5*l4[faulty_id]
    min_error = []
    for k in range(len(indlinks_obs)):
        indlinks = list(indlinks_obs)
        del indlinks[k]
        obs = [g1.indlinks[id] for id in indlinks]
        min_e = np.inf
        for i in [100.0, 1000.0]:
            theta = invopt.solver_mis([g1, g2, g3, g4], [l1[obs], l2[obs], l3[obs], l4[obs]], 
                          indlinks, degree, i*np.ones(degree))
            g1, g2, g3, g4, x1, x2, x3, x4 = get_graphs_linkflows(theta, noisy)
            e = np.linalg.norm(matrix([l1[obs],l2[obs],l3[obs],l4[obs]])-matrix([x1[obs],x2[obs],x3[obs],x4[obs]]))
            if e < min_e: min_e = e
        min_error.append(min_e)
    print min_error
    

def main():
    
    indlinks_obs = [(36,37,1), (13,14,1), (17,8,1), (24,17,1), (28,22,1), (14,13,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    #indlinks_obs = [(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    #indlinks_obs = graph.indlinks.keys()
    #indlinks_obs = [indlinks_obs[3*i] for i in range(len(indlinks_obs)/3)]
    #indlinks_obs = []
    #indlinks_obs = [(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
    
    type = 'Hyperbolic'
    #type = 'Polynomial'
    #test1(type)
    #test2(type)
    #test3(indlinks_obs)
    test4(indlinks_obs, False, type)
    #test5(indlinks_obs, (24,40,1), False)
    
if __name__ == '__main__':
    main()