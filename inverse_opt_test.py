'''
Created on Jun 6, 2014

@author: jeromethai
'''

import numpy as np
import ue_solver as ue
import inverse_opt as invopt
from generate_graph import los_angeles
import matplotlib.pyplot as plt
from cvxopt import matrix
from util import add_noise

a, b = 3.5, 3.0
coef = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
degree = len(coef)
graph = los_angeles(coef, 'Polynomial')[0]
    
    
def display_results(error, true_theta, thetas, delaytype, info):
    """
    Parameters
    ----------
    error: error b/w true link flows and est. link flows
    true_theta: true parameters
    thetas: list of best parameters
    delaytype: type of the delay function
    info: optimal regularization parameter or standard deviation
    """
    xdata = np.linspace(0.0, 2.5, num=100)
    if len(thetas) == 1: label, style = 'estimate', 'b'
    if len(thetas) > 1: label, style = None, 'b--'
    for theta in thetas:
        vals = [1+(theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
        plt.plot(xdata, vals, style, label=label)
    if delaytype == 'Polynomial':
        true_vals = [1+(true_theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    if delaytype == 'Hyperbolic':
        a,b = true_theta
        true_vals = [1 - a/b + a/(b-x) for x in xdata]
    plt.plot( xdata, true_vals, 'r', label='true')
    plt.xlabel('Link flow (1000 veh/h)')
    plt.ylabel('Delay')
    if len(thetas) == 1: plt.title(r'Estimated delay function, error={:.0f}%, beta={:.0e}'.format(100.0*error, info))
    if len(thetas) > 1: plt.title(r'Estimated delay function, mean error={:.0f}%, std={:.3f}'.format(100.0*error, info))
    plt.legend()
    plt.show()


def display_results_multiple(error, ture_theta, thetas, delaytype, std):
    xdata = np.linspace(0.0, 2.5, num=100)
    vals = [1+(best_theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    if delaytype == 'Polynomial':
        true_vals = [1+(true_theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    if delaytype == 'Hyperbolic':
        a,b = true_theta
        true_vals = [1 - a/b + a/(b-x) for x in xdata]

    
def experiment(indlinks_obs, delaytype, betas, noise=0.0, display=False, soft=1000.0):
    """find parameters that minimizes the distance between x^obs_true in NOISY case
    and x^obs generated by each candidate function with PARTIAL observation
    
    Parameters
    ----------
    indlinks_obs: indices of the observed links
    delaytype: type of the delay
    betas: list of smoothing parameters
    noise: std of the noise added to the measured link flows, ff delays, OD demands
    display: if True, display results
    soft: weight put on the observation
    """
    if delaytype == 'Polynomial': true_theta = coef
    if delaytype == 'Hyperbolic': true_theta = (a,b)
    g1, g2, g3, g4 = los_angeles(true_theta, delaytype)
    l1, l2, l3, l4 = ue.solver(g1), ue.solver(g2), ue.solver(g3), ue.solver(g4)
    obs = [g1.indlinks[id] for id in indlinks_obs]
    obs = [int(i) for i in list(np.sort(obs))]
    x1,x2,x3,x4 = l1,l2,l3,l4
    if noise > 0.0:
        x1, x2, x3, x4 = add_noise(l1,noise), add_noise(l2,noise), add_noise(l3,noise), add_noise(l4,noise)
        g1, g2, g3, g4 = los_angeles(true_theta, 'Polynomial', noise)
    theta, xs, beta = invopt.main_solver([g1,g2,g3,g4], [x1[obs],x2[obs],x3[obs],x4[obs]], obs, degree, betas, soft)
    u, v = matrix([l1,l2,l3,l4]), matrix(xs)
    error = np.linalg.norm(u-v, 1) / np.linalg.norm(u, 1)
    if display: display_results(error, true_theta, [theta], delaytype, beta)
    return error, theta


def test1(type, noise):
    """Do multiple experiments for multiple sets of observed links
    
    Parameters:
    -----------
    type: type of the delay function
    noise: std of the noise added to the measured link flows, ff delays, OD demands
    """
    ind_obs = {}
    ind_obs[0] = graph.indlinks.keys()
    ind_obs[1] = [(36,37,1), (13,14,1), (17,8,1), (24,17,1), (28,22,1), (14,13,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    ind_obs[2] = [(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    ind_obs[3] = [(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
    betas = [1e-2, 1e0, 1e2, 1e4, 1e6]
    for k in range(len(ind_obs)): experiment(ind_obs[k], type, betas, noise, True)


def test2(delaytype, noise):
    """Do 20 experiments with random noise and plot the results in the same graph
    
    Parameters
    ----------
    delaytype: type of the delay function
    noise: std of the noise added to the measured link flows, ff delays, OD demands
    """
    if delaytype == 'Polynomial': true_theta = coef
    if delaytype == 'Hyperbolic': true_theta = (a,b)
    ind_obs = {}
    ind_obs[0] = graph.indlinks.keys()
    ind_obs[1] = [(36,37,1), (13,14,1), (17,8,1), (24,17,1), (28,22,1), (14,13,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    ind_obs[2] = [(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    ind_obs[3] = [(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
    for k in range(len(ind_obs)):
        errors, thetas, betas = [], [], [1e-2, 1e0, 1e2, 1e4, 1e6]
        for i in range(20):
            e, theta = experiment(ind_obs[k], delaytype, betas, noise)
            errors.append(e)
            thetas.append(theta)
        mean_error, dev_error = np.mean(errors), np.std(errors)
        display_results(mean_error, true_theta, thetas, delaytype, dev_error)


def test3(delaytype, noise):
    """Do multiple experiments with different values of soft"""
    if delaytype == 'Polynomial': true_theta = coef
    if delaytype == 'Hyperbolic': true_theta = (a,b)
    ind_obs = graph.indlinks.keys()
    for j in [1e-2, 1e0, 1e2, 1e4, 1e6]:
        errors, thetas, betas = [], [], [1e-2, 1e0, 1e2, 1e4, 1e6]
        for i in range(20):
            e, theta = experiment(ind_obs, delaytype, betas, noise, soft=j)
            errors.append(e)
            thetas.append(theta)
        mean_error, dev_error = np.mean(errors), np.std(errors)
        display_results(mean_error, true_theta, thetas, delaytype, dev_error)


def plot_errors():
    """Display errorbars from results of esperiments"""
    x = [100.0*i/60.0 for i in range(5)]
    y1 = [0.002, 4.401, 8.276, 15.655, 21.224]
    y2 = [0.057, 3.222, 6.901, 10.436, 14.643]
    y3 = [6.428, 6.825, 8.050, 9.802, 12.071]
    #y4 = [25.934, ]
    y5 = [0.239, 4.795, 9.810, 12.713, 18.132]
    y6 = [6.431, 6.493, 8.232, 11.177, 12.176]
    y7 = [5.546, 6.223, 7.305, 8.664, 10.949]
    yerr1 = [0, 1.390, 3.273, 4.183, 7.339]
    yerr2 = [0, 1.036, 1.982, 3.265, 6.722]
    yerr3 = [0, 0.736, 1.460, 2.643, 2.804]
    #yerr4 = []
    yerr5 = [0, 1.504, 2.814, 3.420, 6.918]
    yerr6 = [0, 0.899, 1.930, 2.387, 2.150]
    yerr7 = [0, 0.788, 1.280, 1.421, 1.483]
    
    plt.errorbar(x, y1, yerr=yerr1, fmt='-o')
    plt.errorbar(x, y2, yerr=yerr2, fmt='-o')
    plt.errorbar(x, y3, yerr=yerr3, fmt='-o')
    plt.xlabel('Relative amount of standard deviation (%)')
    plt.ylabel('l2-error')
    plt.title('Mean error and standard deviation')
    plt.show()
    
    plt.errorbar(x, y5, yerr=yerr5, fmt='-o')
    plt.errorbar(x, y6, yerr=yerr6, fmt='-o')
    plt.errorbar(x, y7, yerr=yerr7, fmt='-o')
    plt.xlabel('Relative amount of standard deviation (%)')
    plt.ylabel('l2-error')
    plt.title('Mean error and standard deviation')
    plt.show()


def main():
    
    type = 'Polynomial'
    #type = 'Hyperbolic'
    #noise = 0.0
    noise = 0.0
    test1(type, noise)
    #test2(type, noise)
    #test3(type, noise)
    #plot_errors()
    
    
if __name__ == '__main__':
    main()