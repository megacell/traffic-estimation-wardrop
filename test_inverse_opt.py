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
                                              indlinks_obs, degree, smooth, None, max_iter)
    
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
                                                 indlinks_obs, degree, smooth, None, max_iter)
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
    
    
def test3(missing, max_iter, smooth, indlinks_obs=None, soft=None):
    graph1, graph2, graph3 = los_angeles(theta_true, 'Polynomial', True)
    linkflows1 = ue.solver(graph1, update=False)
    linkflows2 = ue.solver(graph2, update=False)
    linkflows3 = ue.solver(graph3, update=False)
    
    if not missing:
        theta = invopt.solver([graph1, graph2, graph3], [linkflows1, linkflows2, linkflows3], degree, smooth)
    else:
        obs = [graph1.indlinks[id] for id in indlinks_obs]
        theta = invopt.solver_mis([graph1, graph2, graph3], [linkflows1[obs], linkflows2[obs], linkflows3[obs]], 
                                  indlinks_obs, degree, smooth, soft, max_iter)
        
    graph1, graph2, graph3 = los_angeles(theta, 'Polynomial', True)
    l_lkflows = [ue.solver(graph1, update=False), ue.solver(graph2, update=False), 
                     ue.solver(graph3, update=False)]
    print 'Estimated parameters'
    print theta
    xdata = np.linspace(0.0, 2.5, num=100)
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
    error = np.linalg.norm(matrix([e1,e2,e3]))
    print error
    #reverse = {v:id for id,v in graph1.indlinks.items()}
    #print [(reverse[i], x1[i], linkflows1[i]) for i in range(n) if abs(e1[i]) > 1e-1]
    #print [(reverse[i], x2[i], linkflows2[i]) for i in range(n) if abs(e2[i]) > 1e-1]
    #print [(reverse[i], x3[i], linkflows3[i]) for i in range(n) if abs(e3[i]) > 1e-1]
    
    plt.plot(xdata, vals, 'r', label='estimate')
    plt.plot( xdata, true_vals, 'b', label='true')
    plt.xlabel('Link flow')
    plt.ylabel('Delay')
    plt.title(r'Estimated delay function. l2-norm error: {:.3f}'.format(error))
    plt.legend()
    plt.show()
    
    
def test4(fvalue=False):
    """3-fold test to find smooth that minimizes the distance between x_true (exact)
    and x generated by each candidate function with TOTAL observation
    results of the 1-fold: best smooth = [600.0, 600.0, 600.0, 600.0, 600.0, 600.0]
    results of the 2-fold: best smooth = [30.0, 30.0, 30.0, 1000.0, 1000.0, 1000.0]
    results of the 3-fold: best smooth = [30.0, 30.0, 60.0, 60.0, 3000.0, 3000.0]
    """
    graph1, graph2, graph3 = los_angeles(theta_true, 'Polynomial', True)
    l1 = ue.solver(graph1, update=False)
    l2 = ue.solver(graph2, update=False)
    l3 = ue.solver(graph3, update=False)
    
    min_error = np.inf
    for i in [10.0, 30.0, 60.0]:
        for j in [30.0, 60.0, 100.0]:
            for k in [1000.0, 3000.0, 6000.0]:
                smooth = np.hstack((i*np.ones(degree/3), j*np.ones(degree/3), k*np.ones(degree/3)))
                if fvalue:
                    theta, f = invopt.solver([graph1, graph2, graph3], [l1, l2, l3], degree, smooth, fvalue=True)
                else:
                    theta = invopt.solver([graph1, graph2, graph3], [l1, l2, l3], degree, smooth)
                g1, g2, g3 = los_angeles(theta, 'Polynomial', True)
                x1 = ue.solver(g1, update=False)
                x2 = ue.solver(g2, update=False) 
                x3 = ue.solver(g3, update=False)
                error = np.linalg.norm(matrix([l1,l2,l3])-matrix([x1,x2,x3]))
                if error < min_error:
                    best_smooth, min_error = (i,j,k), error
                    if fvalue: obj = f
    print best_smooth, min_error
    if fvalue: print obj
    

def test5():
    """3-fold test to find the best smooth using cross validation
    average over the three results with TOTAL observation
    results of the 1-fold best smooth = smooth = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
    results of the 2-fold best smooth = smooth = [30.0, 30.0, 30.0, 1533.3, 1533.3, 1533.3]
    results of the 3-fold best smooth = smooth = [30.0, 30.0, 73.3, 73.3, 3000.0, 3000.0]
    !!!!! For the 3-fold experiment, had to throw away the value of i for the first rotation !!!!!
    """
    graph1, graph2, graph3 = los_angeles(theta_true, 'Polynomial', True)
    l1 = ue.solver(graph1, update=False)
    l2 = ue.solver(graph2, update=False)
    l3 = ue.solver(graph3, update=False)
    
    
    min_error = np.inf
    for i in [10.0, 30.0, 60.0]:
        for j in [30.0, 60.0, 100.0]:
            for k in [1000.0, 3000.0, 6000.0]:
                smooth = np.hstack((i*np.ones(degree/3), j*np.ones(degree/3), k*np.ones(degree/3)))
                theta = invopt.solver([graph1, graph2], [l1, l2], degree, smooth)
                g1, g2, g3 = los_angeles(theta, 'Polynomial', True)
                x3 = ue.solver(g3, update=False)
                error = np.linalg.norm(l3-x3)
                if error < min_error: best_smooth1, min_error = (i,j,k), error
        
    min_error = np.inf
    for i in [10.0, 30.0, 60.0]:
        for j in [30.0, 60.0, 100.0]:
            for k in [1000.0, 3000.0, 6000.0]:
                smooth = np.hstack((i*np.ones(degree/3), j*np.ones(degree/3), k*np.ones(degree/3)))
                theta = invopt.solver([graph1, graph3], [l1, l3], degree, smooth)
                g1, g2, g3 = los_angeles(theta, 'Polynomial', True)
                x2 = ue.solver(g2, update=False)
                error = np.linalg.norm(l2-x2)
                if error < min_error: best_smooth2, min_error = (i,j,k), error
        
    min_error = np.inf
    for i in [10.0, 30.0, 60.0]:
        for j in [60.0, 100.0, 300.0]:
            for k in [1000.0, 3000.0, 6000.0]:
                smooth = np.hstack((i*np.ones(degree/3), j*np.ones(degree/3), k*np.ones(degree/3)))
                theta = invopt.solver([graph3, graph2], [l3, l2], degree, smooth)
                g1, g2, g3 = los_angeles(theta, 'Polynomial', True)
                x1 = ue.solver(g1, update=False)
                error = np.linalg.norm(l1-x1)
                if error < min_error: best_smooth3, min_error = (i,j,k), error
        
    print best_smooth1, best_smooth2, best_smooth3

    
def test6(indlinks_obs, max_iter, fvalue=False):
    """2-fold test to find smooth that minimizes the distance between x^obs_true (exact)
    and x^obs generated by each candidate function with PARTIAL observation
    results for 2-fold with obs = [(17,24),(24,40),(14,21),(16,23)]: best smooth = [30.0, 30.0, 30.0, 1000.0, 1000.0, 1000.0]
    results for 2-fold with obs = [(10,9),(19,18),(4,5),(29,21)]: best smooth = [100.0, 100.0, 100.0, 3000.0, 3000.0, 3000.0]
    """
    graph1, graph2, graph3 = los_angeles(theta_true, 'Polynomial', True)
    l1 = ue.solver(graph1, update=False)
    l2 = ue.solver(graph2, update=False)
    l3 = ue.solver(graph3, update=False)
    obs = [graph1.indlinks[id] for id in indlinks_obs]
    min_error = np.inf
    for i in [30.0, 60.0, 100.0]:
        for j in [600.0, 1000.0, 3000.0]:
            smooth = np.hstack((i*np.ones(degree/2), j*np.ones(degree/2)))
            if fvalue:
                theta, f1, f2 = invopt.solver_mis([graph1, graph2, graph3], [l1[obs], l2[obs], l3[obs]], 
                                  indlinks_obs, degree, smooth, 1000.0, max_iter, fvalue=True)
            else:
                theta = invopt.solver_mis([graph1, graph2, graph3], [l1[obs], l2[obs], l3[obs]], 
                                  indlinks_obs, degree, smooth, 1000.0, max_iter, alt=True)
            g1, g2, g3 = los_angeles(theta, 'Polynomial', True)
            x1 = ue.solver(g1, update=False)
            x2 = ue.solver(g2, update=False) 
            x3 = ue.solver(g3, update=False)
            e = np.linalg.norm(matrix([l1[obs],l2[obs],l3[obs]])-matrix([x1[obs],x2[obs],x3[obs]]))
            if e < min_error:
                best_smooth, min_error, best_theta = (i,j), e, theta
                y1, y2, y3 = x1, x2, x3
                if fvalue: obj1, obj2 = f1, f2
                
    print best_smooth, min_error
    print 'Estimated parameters'
    print best_theta
    if fvalue: print obj1, obj2
    error = np.linalg.norm(matrix([l1, l2, l3])-matrix([y1, y2, y3]))
    xdata = np.linspace(0.0, 2.5, num=100)
    vals = [1+(best_theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    true_vals = [1+(theta_true.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    plt.plot(xdata, vals, 'r', label='estimate')
    plt.plot( xdata, true_vals, 'b', label='true')
    plt.xlabel('Link flow')
    plt.ylabel('Delay')
    plt.title(r'Estimated delay function. l2-norm error: {:.3f}'.format(error))
    plt.legend()
    plt.show()
    
    
def test7(indlinks_obs, max_iter, fvalue=False):
    """2-fold test to find smooth that minimizes the distance between x^obs_true (exact)
    and x^obs generated by each candidate function in NOISY case
    results for 2-fold with all obs: best smooth = [30.0, 30.0, 30.0, 1000.0, 1000.0, 1000.0]
    results for 2-fold with obs = [(17,24),(24,40),(14,21),(16,23)]: best smooth = [30.0, 30.0, 30.0, 1000.0, 1000.0, 1000.0]
    results for 2-fold with obs = [(10,9),(19,18),(4,5),(29,21)]: best smooth = [10.0, 10.0, 10.0, 3000.0, 3000.0, 3000.0]
    """
    graph1, graph2, graph3 = los_angeles(theta_true, 'Polynomial', True)
    l1 = ue.solver(graph1, update=False)
    l2 = ue.solver(graph2, update=False)
    l3 = ue.solver(graph3, update=False)
    obs = [graph1.indlinks[id] for id in indlinks_obs]
    
    graph1, graph2, graph3 = los_angeles(theta_true, 'Polynomial', True, True)
    z1, z2, z3 = l1, l2, l3
    np.random.seed(21)
    l1, l2, l3 = matrix(np.random.normal(l1, l1/30.0)), matrix(np.random.normal(l2, l2/30.0)), matrix(np.random.normal(l3, l3/30.0))
    
    min_error = np.inf
    for i in [30.0, 60.0, 100.0]:
        for j in [600.0, 1000.0, 3000.0, 6000.0]:
            smooth = np.hstack((i*np.ones(degree/2), j*np.ones(degree/2)))
            if fvalue:
                theta, f1, f2 = invopt.solver_mis([graph1, graph2, graph3], [l1[obs], l2[obs], l3[obs]], 
                                  indlinks_obs, degree, smooth, 1000.0, max_iter, fvalue=True)
            else:
                theta = invopt.solver_mis([graph1, graph2, graph3], [l1[obs], l2[obs], l3[obs]], 
                                  indlinks_obs, degree, smooth, 1000.0, max_iter, alt=True)
            g1, g2, g3 = los_angeles(theta, 'Polynomial', True, True)
            x1 = ue.solver(g1, update=False)
            x2 = ue.solver(g2, update=False) 
            x3 = ue.solver(g3, update=False)
            e = np.linalg.norm(matrix([l1[obs],l2[obs],l3[obs]])-matrix([x1[obs],x2[obs],x3[obs]]))
            if e < min_error:
                best_smooth, min_error, best_theta = (i,j), e, theta
                y1, y2, y3 = x1, x2, x3
                if fvalue: obj1, obj2 = f1, f2
                
    print best_smooth, min_error
    print 'Estimated parameters'
    print best_theta
    if fvalue: print obj1, obj2
    l1, l2, l3 = z1, z2, z3
    error = np.linalg.norm(matrix([l1, l2, l3])-matrix([y1, y2, y3]))
    xdata = np.linspace(0.0, 2.5, num=100)
    vals = [1+(best_theta.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    true_vals = [1+(theta_true.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    plt.plot(xdata, vals, 'r', label='estimate')
    plt.plot( xdata, true_vals, 'b', label='true')
    plt.xlabel('Link flow')
    plt.ylabel('Delay')
    plt.title(r'Estimated delay function. l2-norm error: {:.3f}'.format(error))
    plt.legend()
    plt.show()
    
    
def test8(max_iter):
    """Cross validation to find a sensor that has been attacked
    results when obs = [(17,24),(24,40),(14,21),(16,23)] and (24,40) has been attacked
    """
    indlinks_obs = [(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    graph1, graph2, graph3 = los_angeles(theta_true, 'Polynomial', multiple=True)
    l1 = ue.solver(graph1, update=False)
    l2 = ue.solver(graph2, update=False)
    l3 = ue.solver(graph3, update=False)
    faulty_id = graph1.indlinks[(24,40,1)]
    l1[faulty_id] = 0.7*l1[faulty_id]
    l2[faulty_id] = 0.7*l2[faulty_id]
    l3[faulty_id] = 0.7*l3[faulty_id]
    
    min_error = []
    for k in range(4):
        indlinks = list(indlinks_obs)
        del indlinks[k]
        obs = [graph1.indlinks[id] for id in indlinks]
        min_e = np.inf
        for i in [10.0, 30.0, 60.0]:
            for j in [600.0, 1000.0, 3000.0]:
                smooth = np.hstack((i*np.ones(degree/2), j*np.ones(degree/2)))
                theta = invopt.solver_mis([graph1, graph2, graph3], [l1[obs], l2[obs], l3[obs]], 
                                  indlinks, degree, smooth, None, max_iter)
                g1, g2, g3 = los_angeles(theta, 'Polynomial', True)
                x1 = ue.solver(g1, update=False)
                x2 = ue.solver(g2, update=False) 
                x3 = ue.solver(g3, update=False)
                e = np.linalg.norm(matrix([l1[obs],l2[obs],l3[obs]])-matrix([x1[obs],x2[obs],x3[obs]]))
                if e < min_e: min_e = e
        min_error.append(min_e)
    print min_error
        
    
def test9(max_iter):
    """Cross validation to find a sensor that has been attacked
    results when obs = [(17,24),(24,40),(14,21),(16,23)] and (24,40) has been attacked
    in the NOISY case
    """
    indlinks_obs = [(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    graph1, graph2, graph3 = los_angeles(theta_true, 'Polynomial', multiple=True)
    l1 = ue.solver(graph1, update=False)
    l2 = ue.solver(graph2, update=False)
    l3 = ue.solver(graph3, update=False)
    np.random.seed(21)
    l1, l2, l3 = matrix(np.random.normal(l1, l1/30.0)), matrix(np.random.normal(l2, l2/30.0)), matrix(np.random.normal(l3, l3/30.0))
    faulty_id = graph1.indlinks[(24,40,1)]
    l1[faulty_id] = 0.3*l1[faulty_id]
    l2[faulty_id] = 0.3*l2[faulty_id]
    l3[faulty_id] = 0.3*l3[faulty_id]
    graph1, graph2, graph3 = los_angeles(theta_true, 'Polynomial', multiple=True, noisy=True)
    min_error = []
    for k in range(4):
        indlinks = list(indlinks_obs)
        del indlinks[k]
        obs = [graph1.indlinks[id] for id in indlinks]
        min_e = np.inf
        for i in [10.0, 30.0, 60.0]:
            for j in [600.0, 1000.0, 3000.0]:
                smooth = np.hstack((i*np.ones(degree/2), j*np.ones(degree/2)))
                theta = invopt.solver_mis([graph1, graph2, graph3], [l1[obs], l2[obs], l3[obs]], 
                                  indlinks, degree, smooth, 1000.0, max_iter)
                g1, g2, g3 = los_angeles(theta, 'Polynomial', multiple=True, noisy=True)
                x1 = ue.solver(g1, update=False)
                x2 = ue.solver(g2, update=False) 
                x3 = ue.solver(g3, update=False)
                e = np.linalg.norm(matrix([l1[obs],l2[obs],l3[obs]])-matrix([x1[obs],x2[obs],x3[obs]]))
                if e < min_e: min_e = e
        min_error.append(min_e)
    print min_error
    

def main():
    
    #indlinks_obs = [(8,17,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1), (23,28,1), (12,34,1), (13,14,1)]
    indlinks_obs = [(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    #graph = los_angeles(theta_true, 'Polynomial')
    #indlinks_obs = graph.indlinks.keys()
    #indlinks_obs = []
    #indlinks_obs = [(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
    
    #smooth = 500.0*np.ones(degree)
    #smooth = np.hstack((100.0*np.ones(degree/2), 3000.0*np.ones(degree/2)))
    #smooth = np.hstack((30.0*np.ones(degree/3), 60.0*np.ones(degree/3), 3000.0*np.ones(degree/3)))
    
    #test3(True, 10, smooth, indlinks_obs, 0.0)
    #test4(True)
    #test5()
    test6(indlinks_obs, 10, True)
    #test7(indlinks_obs, 10, False)
    #test9(10)
    
if __name__ == '__main__':
    main()