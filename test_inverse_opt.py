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
theta_true = matrix([0.0, 1.0, 0.0])
theta_true /= np.sum(theta_true)
degree = len(theta_true)

def test1(missing):
    
    graph1 = small_grid(od_flows1, 'Polynomial', theta_true)
    graph2 = small_grid(od_flows2, 'Polynomial', theta_true)
    linkflows1 = ue.solver(graph1, update=False)
    linkflows2 = ue.solver(graph2, update=False)
    
    if not missing:
        
        theta = invopt.solver([graph1, graph2], [linkflows1, linkflows2], degree)
    
    else:
    
        indlinks_obs = [(2,1,1), (6,5,1), (3,2,1), (5,4,1)]
        obs = [graph1.indlinks[id] for id in indlinks_obs]
        theta = invopt.solver_mis([graph1, graph2], [linkflows1[obs], linkflows2[obs]], indlinks_obs, degree, 2)
    
    print 'Estimated parameters'
    print theta
    xdata = np.linspace(0.0, 5.0, num=10)
    vals = [1+(theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    true_vals = [1+(theta_true.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    scale = sum(true_vals) / sum(vals)
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
    
    for k in range(1, max_obs+1):
        
        tmp2 = 0
        for i in range(trials):
            ind = range(n)
            obs = np.random.permutation(n)[range(k)]
            indlinks_obs = [indlinks2[i] for i in obs]
            theta = invopt.solver_mis([graph1, graph2], [linkflows1[matrix(obs)], linkflows2[matrix(obs)]], indlinks_obs, degree, 2)
            x1 = ue.solver(small_grid(od_flows1, 'Polynomial', theta), update=False)
            x2 = ue.solver(small_grid(od_flows2, 'Polynomial', theta), update=False)
            tmp1 = matrix([linkflows1, linkflows2]) - matrix([x1, x2])
            tmp2 += tmp1.T*tmp1
        tmp2 /= float(trials)
        error.append(tmp2[0])
    
    plt.plot(range(1,max_obs+1), error)
    plt.xlabel('Number of sensors')
    plt.ylabel('Error')
    plt.title(r'Error in estimate')
    plt.legend()
    plt.show()
    
def main():
    #test1(True)
    test2(8, 20)
    
    
if __name__ == '__main__':
    main()