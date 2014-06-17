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
theta_true = matrix([0.0, 0.0, 0.0, 1.0, 0.0, 0.0])
theta_true /= np.sum(theta_true)

theta_true *= 0.15
degree = len(theta_true)

def test1(missing, max_iter, alpha=None):
    
    graph1 = small_grid(od_flows1, 'Polynomial', theta_true)
    graph2 = small_grid(od_flows2, 'Polynomial', theta_true)
    linkflows1 = ue.solver(graph1, update=False)
    linkflows2 = ue.solver(graph2, update=False)
    
    if not missing:
        
        theta = invopt.solver([graph1, graph2], [linkflows1, linkflows2], degree, smooth, alpha)
    
    else:
    
        #indlinks_obs = [(6,5,1)]
        #indlinks_obs = []
        #indlinks_obs = [(2, 1, 1), (5, 4, 1), (1, 4, 1), (2, 5, 1), (3, 2, 1), (3, 6, 1), (5, 2, 1), (6, 5, 1)]
        #indlinks_obs = [(2, 1, 1), (5, 4, 1), (1, 4, 1), (2, 5, 1), (5, 2, 1), (6, 5, 1)]
        indlinks_obs = [(3, 2, 1), (3, 6, 1), (2, 1, 1), (2, 5, 1)]
        #indlinks_obs = graph1.indlinks.keys()
        obs = [graph1.indlinks[id] for id in indlinks_obs]
        theta, l_lkflows = invopt.solver_mis([graph1, graph2], [linkflows1[obs], linkflows2[obs]],
                                              indlinks_obs, degree, smooth, alpha, max_iter)
    
    graph1 = small_grid(od_flows1, 'Polynomial', theta)
    graph2 = small_grid(od_flows2, 'Polynomial', theta)
    l_lkflows = [ue.solver(graph1, update=False), ue.solver(graph2, update=False)]
    print 'Estimated parameters'
    print theta
    xdata = np.linspace(0.0, 5.0, num=10)
    vals = [1+(theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    true_vals = [1+(theta_true.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    scale = sum(true_vals) / sum(vals)
    #scale = 1
    scaled_vals = [scale*val for val in vals]
    
    print (l_lkflows[0]-linkflows1).T*(l_lkflows[0]-linkflows1)
    print (l_lkflows[1]-linkflows2).T*(l_lkflows[1]-linkflows2)
    
    plt.plot(xdata, scaled_vals, 'r', label='estimate')
    plt.plot( xdata, true_vals, 'b', label='true')
    plt.xlabel('Link flow')
    plt.ylabel('Delay')
    plt.title(r'Estimated delay function')
    plt.legend()
    plt.show()
    
    
def test2(max_obs, trials, max_iter):
    
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
            theta, l_lkflows = invopt.solver_mis([graph1, graph2], [linkflows1[matrix(obs)], linkflows2[matrix(obs)]], 
                                                 indlinks_obs, degree, smooth, alpha, max_iter=max_iter)
            tmp1 = matrix([linkflows1[mis], linkflows2[mis]]) - matrix([l_lkflows[0][mis], l_lkflows[1][mis]])
            tmp2.append(max(abs(tmp1)) / float(n-k))
            
        error.append(np.mean(tmp2))
        error_var.append(np.std(tmp2))
    
    plt.plot(range(max_obs+1), error)
    plt.errorbar(range(max_obs+1), error, yerr = error_var, fmt = 'o')
    plt.xlabel('Number of sensors')
    plt.ylabel('Error')
    plt.title(r'Error in estimate')
    plt.legend()
    plt.show()
    
    
def test3(missing, max_iter, smooth):
    graph1, graph2, graph3 = los_angeles(theta_true, 'Polynomial', True)
    linkflows1 = ue.solver(graph1)
    linkflows2 = ue.solver(graph2)
    linkflows3 = ue.solver(graph3)
    
    if not missing:
        theta = invopt.solver([graph1, graph2, graph3], [linkflows1, linkflows2, linkflows3], degree, smooth)
    else:
        indlinks_obs = [(8,17,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1), (23,28,1), (12,34,1), (13,14,1)]
        #indlinks_obs = graph1.indlinks.keys()
        #indlinks_obs = []
        #indlinks_obs = [(41,8,1), (42,17,1), (43,24,1), (5,12,1), (35,36,1), (1,2,1)]
        obs = [graph1.indlinks[id] for id in indlinks_obs]
        theta, l_lkflows = invopt.solver_mis([graph1, graph2, graph3], [linkflows1[obs], linkflows2[obs], linkflows3[obs]], 
                                  indlinks_obs, degree, smooth, max_iter)
        
    graph1, graph2, graph3 = los_angeles(theta, 'Polynomial', True)
    l_lkflows = [ue.solver(graph1, update=False), ue.solver(graph2, update=False), 
                     ue.solver(graph3, update=False)]
    print 'Estimated parameters'
    print theta
    xdata = np.linspace(0.0, 3.0, num=10)
    vals = [1+(theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    true_vals = [1+(theta_true.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    
    n = graph1.numlinks
    e1 = abs(linkflows1-l_lkflows[0])
    e2 = abs(linkflows2-l_lkflows[1])
    e3 = abs(linkflows3-l_lkflows[2])
    #print l_lkflows[2][graph1.indlinks[(15,22,1)]]
    print 'l1 norm'
    print np.mean(e1)
    print np.mean(e2)
    print np.mean(e3)
    print 'l2 norm'
    print np.linalg.norm(e1)
    print np.linalg.norm(e2)
    print np.linalg.norm(e3)
    #reverse = {v:id for id,v in graph1.indlinks.items()}
    #print [(reverse[i], x1[i], linkflows1[i]) for i in range(n) if abs(e1[i]) > 1e-1]
    #print [(reverse[i], x2[i], linkflows2[i]) for i in range(n) if abs(e2[i]) > 1e-1]
    #print [(reverse[i], x3[i], linkflows3[i]) for i in range(n) if abs(e3[i]) > 1e-1]
    
    
    plt.plot(xdata, vals, 'r', label='estimate')
    plt.plot( xdata, true_vals, 'b', label='true')
    plt.xlabel('Link flow')
    plt.ylabel('Delay')
    plt.title(r'Estimated delay function')
    plt.legend()
    plt.show()
    
    
def main():
    #test1(True, 5, 1.5*alpha)
    #test2(7, 20, 3)
    test3(True, 2, 600.0*np.ones(degree))
    
    
if __name__ == '__main__':
    main()