'''
Created on Jun 6, 2014

@author: jeromethai
'''

import numpy as np
import ue_solver as ue
import inverse_opt as invopt
from test_graph import small_grid
import matplotlib.pyplot as plt
from cvxopt import matrix

od_flows1 = [3.0, 3.0, 1.0, 1.0];
od_flows2 = [1.0, 1.0, 1.0, 4.0];
theta_true = matrix([3.0, 2.0, 1.0])
theta_true /= np.sum(theta_true)
degree = len(theta_true)

def main(missing):
    
    graph1 = small_grid(od_flows1, 'Polynomial', theta_true)
    graph2 = small_grid(od_flows2, 'Polynomial', theta_true)
    linkflows1 = ue.solver(graph1, update=False)
    linkflows2 = ue.solver(graph2, update=False)
    
    if not missing:
        
        theta = invopt.solver([graph1, graph2], [linkflows1, linkflows2], degree)
        print 'Estimated parameters'
        print theta
    
    else:
    
        indlinks_obs, obs = [], []
        for id in indlinks_obs: obs.append(graph1.indlinks[id])
        theta = invopt.solver_mis([graph1, graph2], [linkflows1[obs], linkflows2[obs]], indlinks_obs, degree, 1)
    
    xdata = np.linspace(0.0, 5.0, num=10)
    vals = [1+(theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    true_vals = [1+(theta_true.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    scale = sum(true_vals) / sum(vals)
    scaled_vals = [scale*val for val in vals]
    
    plt.plot(xdata, scaled_vals, 'r', label='estimate')
    plt.plot( xdata, true_vals, 'b', label='true')
    plt.xlabel('Link flow')
    plt.ylabel('Delay')
    plt.title(r'Estimated delay function')
    plt.legend()
    plt.show()
    
    
if __name__ == '__main__':
    main(True)