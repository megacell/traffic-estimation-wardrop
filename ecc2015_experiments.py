'''
Created on Oct 3, 2014

@author: jeromethai

Experiments for ECC 2015
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


def experiment(indlinks_obs, delaytype, display=False):
    """Generate the graph of L.A.
    
    Parameters
    ----------
    indlinks_obs = indices of the observed links
    delaytype: type of the delay
    display: if True, display results
    """
    if delaytype == 'Polynomial': true_theta = coef
    if delaytype == 'Hyperbolic': true_theta = (a,b)
    g1, g2, g3, g4 = los_angeles(true_theta, delaytype)
    x1, x2, x3, x4 = ue.solver(g1), ue.solver(g2), ue.solver(g3), ue.solver(g4)
    obs = [g1.indlinks[id] for id in indlinks_obs]
    obs = [int(i) for i in list(np.sort(obs))]
    w_multi = [.01, .1, .5, .9, .99]
    w_reg = [1e-6, 1e-4, 1e-2, 1e0, 1e2]
    r_gap, r_obs, x_est = invopt.multi_objective_solver([g1,g2,g3,g4], [x1[obs],x2[obs],x3[obs],x4[obs]], obs, degree, w_multi, w_reg)
    print r_gap
    print r_obs
    m,n = r_gap.size
    print m
    print n
    print x_est
    r_est = matrix(0.0, (m,n))
    u = matrix([x1,x2,x3,x4])
    for i in range(m):
        for j in range(n):
            r_est[i,j] = np.linalg.norm(u-x_est[(i,j)], 1) / np.linalg.norm(u, 1)
    print r_est
    
    
def display_results():
    """"""
    
    
    

def test1(type):
    ind_obs = {}
    ind_obs[0] = [(36,37,1), (13,14,1), (17,8,1), (24,17,1), (28,22,1), (14,13,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    #ind_obs[2] = [(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    #ind_obs[3] = [(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
    for k in range(len(ind_obs)): experiment(ind_obs[k], type)


def main():
    type = 'Hyperbolic'
    test1(type)


if __name__ == '__main__':
    main()