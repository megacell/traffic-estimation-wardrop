'''
Created on Sep 9, 2014

@author: jeromethai
'''

import inverse_opt as invopt
import ue_solver as ue
import toll_pricing as tp
from generate_graph import los_angeles
from cvxopt import matrix
from inverse_opt_test import experiment


def test_toll_pricing(theta, display=False):
    """Test the toll pricing model
    1. compute optimal toll
    2. compute tolled ue, ue, so total travel times
    3. display results
    """
    g1, g2, g3, g4 = los_angeles(theta, 'Polynomial')
    tolls, tolls_collected, costs, ue_costs, so_costs, weights = [], [], [], [], [], []
    for graph in [g2,g3,g4]:
        ffdelays, slopes = graph.get_ffdelays(), graph.get_slopes()
        Aeq, beq = ue.constraints(graph)
        coefs = invopt.compute_coefs(ffdelays, slopes, theta)
        data = (Aeq, beq, ffdelays, coefs, 'Polynomial')
        toll, cost, toll_collected, weight = tp.main_solver(graph, theta)
        tolls.append(toll)
        tolls_collected.append(toll_collected)
        weights.append(weight)
        costs.append(tp.compute_cost(data, toll)[0]) # tolled equilibrium
        ue_costs.append(tp.compute_cost(data)[0]) # ue equilibrium
        so_costs.append(tp.compute_cost(data, 0.0, True)[0]) # so equilibrium
    if display:
        print 'toll collected:', tolls_collected
        print 'SO total travel time:', so_costs
        print 'tolled total travel time:', costs
        print 'UE total travel time:', ue_costs
        print 'weight:', weights
    return tolls, tolls_collected, so_costs, costs, ue_costs


def two_step_test():
    """"""
    graph = los_angeles()[3]
    ind_obs = {}
    ind_obs[0] = graph.indlinks.keys()
    #ind_obs[1] = [(36,37,1), (13,14,1), (17,8,1), (24,17,1), (28,22,1), (14,13,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    #ind_obs[2] = [(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    #ind_obs[3] = [(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
    betas = [1e-2, 1e0, 1e2, 1e4, 1e6]
    for k in range(len(ind_obs)):
        error, theta = experiment(ind_obs[k], 'Hyperbolic')
        tolls, tolls_collected, so_costs, costs, ue_costs = test_toll_pricing(theta, True)
        #tolls_collected2, so_costs2, costs2, ue_costs2 = [], [], [], []
        #for toll, data in zip(tolls, datas):
        #    data = (data[0], data[1], data[2], data[3], 'Hyperbolic')
        #    costs.append(tp.compute_cost(data, toll)[0]) # tolled equilibrium
        #    ue_costs.append(tp.compute_cost(data)[0]) # ue equilibrium
        #    so_costs.append(tp.compute_cost(data, 0.0, True)[0]) # so equilibrium


def main():
    theta = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
    #test_toll_pricing(theta, True)
    two_step_test()


if __name__ == '__main__':
    main()