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
theta_true = matrix([1.0, 2.0, 3.0])
theta_true /= np.sum(theta_true)
degree = 3

def main():
    graph1 = small_grid(od_flows1, 'Polynomial', theta_true)
    graph2 = small_grid(od_flows2, 'Polynomial', theta_true)
    linkflows1 = ue.solver(graph1, update=False)
    linkflows2 = ue.solver(graph2, update=False)
    c, A, b, Aeq, beq = invopt.constraints([graph1, graph2], [linkflows1, linkflows2], degree)
    theta = invopt.solver([graph1, graph2], [linkflows1, linkflows2], degree)
    print 'Estimated parameters'
    print theta
    
    xdata = np.linspace(0.0, 5.0, num=10)
    vals = [(theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    true_vals = [(theta_true.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    scale = sum(true_vals) / sum(vals)
    #scale = 1
    scaled_vals = [scale*val for val in vals]
    
    #Check if we find the right link flows
    x1 = ue.solver(small_grid(od_flows1, 'Polynomial', theta), update=False)
    x2 = ue.solver(small_grid(od_flows2, 'Polynomial', theta), update=False)
    print linkflows1
    print x1
    print linkflows2
    print x2
    
    plt.plot(xdata, scaled_vals, 'r', label='estimate')
    plt.plot( xdata, true_vals, 'b', label='true')
    plt.xlabel('Link flow')
    plt.ylabel('Delay')
    plt.title(r'Estimated delay function')
    plt.legend()
    plt.show()
    
    
if __name__ == '__main__':
    main()