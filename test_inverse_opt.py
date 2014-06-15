'''
Created on Jun 6, 2014

@author: jeromethai
'''

import numpy as np
import ue_solver as ue
import inverse_opt as invopt
from test_graph import small_grid, los_angeles
import matplotlib.pyplot as plt
from cvxopt import matrix
import draw_graph as d

od_flows1 = [3.0, 3.0, 1.0, 1.0];
od_flows2 = [1.0, 1.0, 1.0, 4.0];
theta_true = matrix([0.0, 0.0, 1.0, 1.0])
theta_true /= np.sum(theta_true)

alpha = 0.15
theta_true *= alpha
degree = len(theta_true)

def test1(missing):
    
    graph1 = small_grid(od_flows1, 'Polynomial', theta_true)
    graph2 = small_grid(od_flows2, 'Polynomial', theta_true)
    linkflows1 = ue.solver(graph1, update=False)
    linkflows2 = ue.solver(graph2, update=False)
    
    if not missing:
        
        theta = invopt.solver([graph1, graph2], [linkflows1, linkflows2], degree)
    
    else:
    
        indlinks_obs = [(2,1,1), (6,5,1), (3,2,1), (5,4,1), (2,5,1), (1,4,1), (3,6,1)]
        #indlinks_obs = []
        obs = [graph1.indlinks[id] for id in indlinks_obs]
        theta = invopt.solver_mis([graph1, graph2], [linkflows1[obs], linkflows2[obs]], indlinks_obs, degree, max_iter=5)
    
    print 'Estimated parameters'
    print theta
    xdata = np.linspace(0.0, 5.0, num=10)
    vals = [1+(theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    true_vals = [1+(theta_true.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    #scale = sum(true_vals) / sum(vals)
    scale = 1
    scaled_vals = [scale*val for val in vals]
    
    x1 = ue.solver(small_grid(od_flows1, 'Polynomial', theta), update=False)
    x2 = ue.solver(small_grid(od_flows2, 'Polynomial', theta), update=False)
    print (linkflows1-x1).T*(linkflows1-x1)
    print (linkflows2-x2).T*(linkflows2-x2)
    
    plt.plot(xdata, scaled_vals, 'r', label='estimate')
    plt.plot( xdata, true_vals, 'b', label='true')
    plt.xlabel('Link flow')
    plt.ylabel('Delay')
    plt.title(r'Estimated delay function')
    plt.legend()
    plt.show()
    
    
def test2(max_obs, trials):
    
    graph1 = small_grid(od_flows1, 'Polynomial', theta_true)
    graph2 = small_grid(od_flows2, 'Polynomial', theta_true)
    linkflows1 = ue.solver(graph1, update=False)
    linkflows2 = ue.solver(graph2, update=False)
    n = graph1.numlinks
    indlinks = graph1.indlinks
    indlinks2 = {value: key for key, value in indlinks.items()}
    
    error = []
    error_var = []
    
    for k in range(max_obs+1):
        
        tmp2 = []
        
        for i in range(trials):
            obs = np.random.permutation(n)[range(k)]
            mis = range(n)
            for j in obs: mis.remove(j)
            indlinks_obs = [indlinks2[i] for i in obs]
            theta = invopt.solver_mis([graph1, graph2], [linkflows1[matrix(obs)], linkflows2[matrix(obs)]], indlinks_obs, degree, max_iter=5)
            x1 = ue.solver(small_grid(od_flows1, 'Polynomial', theta), update=False)
            x2 = ue.solver(small_grid(od_flows2, 'Polynomial', theta), update=False)
            tmp1 = matrix([linkflows1[mis], linkflows2[mis]]) - matrix([x1[mis], x2[mis]])
            tmp2.append(((tmp1.T*tmp1) / float(n-k))[0])
            
        error.append(np.mean(tmp2))
        error_var.append(np.std(tmp2))
    
    plt.plot(range(max_obs+1), error)
    plt.errorbar(range(max_obs+1), error, yerr = error_var, fmt = 'o')
    plt.xlabel('Number of sensors')
    plt.ylabel('Error')
    plt.title(r'Error in estimate')
    plt.legend()
    plt.show()
    
    
def test3(missing):
    graph1, graph2, graph3 = los_angeles(theta_true, 'Polynomial', True)
    linkflows1 = ue.solver(graph1)
    linkflows2 = ue.solver(graph2)
    linkflows3 = ue.solver(graph3)
    
    if not missing:
        theta = invopt.solver([graph1, graph2, graph3], [linkflows1, linkflows2, linkflows3], degree)
    else:
        indlinks_obs = []
        obs = [graph1.indlinks[id] for id in indlinks_obs]
        theta = invopt.solver_mis([graph1, graph2, graph3], [linkflows1[obs], linkflows2[obs], linkflows3[obs]], 
                                  indlinks_obs, degree, max_iter=5)
    
    print 'Estimated parameters'
    print theta
    xdata = np.linspace(0.0, 5.0, num=10)
    vals = [1+(theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    true_vals = [1+(theta_true.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    
    graph1, graph2, graph3 = los_angeles(theta, 'Polynomial', True)
    x1 = ue.solver(graph1, update=False)
    x2 = ue.solver(graph2, update=False)
    x3 = ue.solver(graph3, update=False)
    print (linkflows1-x1).T*(linkflows1-x1)
    print (linkflows2-x2).T*(linkflows2-x2)
    print (linkflows3-x3).T*(linkflows3-x3)
    
    plt.plot(xdata, vals, 'r', label='estimate')
    plt.plot( xdata, true_vals, 'b', label='true')
    plt.xlabel('Link flow')
    plt.ylabel('Delay')
    plt.title(r'Estimated delay function')
    plt.legend()
    plt.show()
    
    
def main():
    #test1(True)
    #test2(7, 50)
    test3(True)
    
    
if __name__ == '__main__':
    main()