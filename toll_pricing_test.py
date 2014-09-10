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


def test_toll_pricing(theta, display=False, ws=[1e-8, 1e-6, 1e-4]):
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
        pm = invopt.compute_coefs(ffdelays, slopes, theta)
        data = (Aeq, beq, ffdelays, pm, 'Polynomial')
        toll, cost, toll_collected, weight = tp.main_solver(graph, theta, ws)
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


def two_step_test(ind_obs):
    """"""
    error, theta = experiment(ind_obs, 'Hyperbolic', display=True)
    tolls, tolls_collected, so_costs, costs, ue_costs = test_toll_pricing(theta, True, [1e-6])
    tolls_collected2, so_costs2, costs2, ue_costs2 = [], [], [], []
    g1, g2, g3, g4 = los_angeles((3.5, 3.0), 'Hyperbolic')
    for toll, graph in zip(tolls, [g2,g3,g4]):
        Aeq, beq, ffdelays, pm, type = ue.get_data(graph)
        data = (Aeq, beq, ffdelays, pm, type)
        cost, x = tp.compute_cost(data, toll)
        costs2.append(cost) # tolled equilibrium
        tolls_collected2.append((toll.T*x)[0])
        ue_costs2.append(tp.compute_cost(data)[0]) # ue equilibrium
        so_costs2.append(tp.compute_cost(data, 0.0, True)[0]) # so equilibrium
    print 'toll collected:', tolls_collected, tolls_collected2
    print 'SO total travel time:', so_costs, so_costs2
    print 'tolled total travel time:', costs, costs2
    print 'UE total travel time:', ue_costs, ue_costs2


def plot_results():
    """Plot results from two_step_test()"""
    # results for ind_obs=graph.indlinks.keys()
    tolls_collected = [[161.15742640370462, 416.80264009965265, 948.124688583883], [161.4013529462645, 417.25218843773354, 950.1674351402539]]
    so_costs = [[2519.5329969919685, 3679.4257321044365, 5325.940489651684], [2544.920874812245, 3718.292867749051, 5383.190501656683]]
    costs = [[2525.7590340322095, 3679.4261220632625, 5356.8071818038225], [2551.0113714206195, 3718.3028736461647, 5412.310857035818]]
    ue_costs = [[2548.8303153671195, 3736.67165067958, 5444.499035875823], [2574.155371287541, 3775.672784740404, 5503.564654849874]]
    # results for ind_obs=[(36,37,1), (13,14,1), (17,8,1), (24,17,1), (28,22,1), (14,13,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    tolls_collected = [[355.56391097372205, 701.3544774585599, 806.6640717434449], [374.34816990991186, 736.3763770517197, 842.0768966826593]]
    so_costs = [[2174.1387842065237, 3306.974167934956, 4969.5560792612605], [2544.920874812245, 3718.292867749051, 5383.190501656683]]
    costs = [[2201.273656230071, 3338.7536674213156, 4988.732654672769], [2564.0153714222406, 3746.115605602468, 5398.046518880896]]
    ue_costs = [[2209.622435372754, 3380.7186533697286, 5104.755602389524], [2574.155371287541, 3775.672784740404, 5503.564654849874]]
    # results for ind_obs=[(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    tolls_collected = [[339.37838703976513, 737.9887571945254, 930.0531221139565], [353.28118085588176, 763.3708450755438, 955.644104623869]]
    so_costs = [[2520.1983680076073, 3833.234908857907, 5645.320052639795], [2544.920874812245, 3718.292867749051, 5383.190501656683]]
    costs = [[2536.138622203524, 3862.614483834193, 5668.430527845267], [2562.897499637454, 3750.9322754235955, 5406.402328441233]]
    ue_costs = [[2556.0266644367803, 3903.9758481975805, 5741.545038028959], [2574.155371287541, 3775.672784740404, 5503.564654849874]]
    # results for ind_obs=[(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
    tolls_collected = [[443.85669649193784, 837.5428841560781, 1380.117295219297], [479.57423528502, 931.2309202734045, 1499.6972922042573]]
    so_costs = [[2158.995768390798, 3843.4202195503713, 7595.092411648515], [2544.920874812245, 3718.292867749051, 5383.190501656683]]
    costs = [[2192.900868595389, 3888.0172745668465, 7626.17585256234], [2575.5748442305207, 3755.844017718783, 5450.304046277645]]
    ue_costs = [[2231.1674179937672, 3999.772800567613, 7867.415688946321], [2574.155371287541, 3775.672784740404, 5503.564654849874]]


def main():
    theta = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
    graph = los_angeles()[3]
    ind_obs = {}
    ind_obs[0] = graph.indlinks.keys()
    ind_obs[1] = [(36,37,1), (13,14,1), (17,8,1), (24,17,1), (28,22,1), (14,13,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    ind_obs[2] = [(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    ind_obs[3] = [(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
    
    #test_toll_pricing(theta, True)
    two_step_test(ind_obs[3])
    


if __name__ == '__main__':
    main()