'''
Created on Sep 9, 2014

@author: jeromethai
'''

import numpy as np
import inverse_opt as invopt
import ue_solver as ue
import toll_pricing as tp
from generate_graph import los_angeles
from cvxopt import matrix
from inverse_opt_test import experiment
import matplotlib.pyplot as plt


def test_toll_pricing(theta, display=False, ws=[1e-8, 1e-6, 1e-4]):
    """Test the toll pricing model
    1. compute optimal toll
    2. compute tolled_ue, ue, so total travel times
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


def two_step_test(ind_obs, delaytype):
    """Estimate then determine the optimal toll pricing
    1. Estimate the delay function 
    2. Compute the optimal toll pricing using the estimated latency function
    3. Compute the real effect of the toll pricing on the true model
    """
    # 1. Estimate
    error, theta = experiment(ind_obs, delaytype, display=True)
    # 2. Compute the optimal toll pricing using the estimated latency function
    tolls, tolls_collected, so_costs, costs, ue_costs = test_toll_pricing(theta, True, [1e-6])
    # 3. Compute the real effect of the toll pricing on the true model
    tolls_collected2, so_costs2, costs2, ue_costs2 = [], [], [], []
    if delaytype == 'Hyperbolic': g1, g2, g3, g4 = los_angeles((3.5, 3.0), 'Hyperbolic')
    if delaytype == 'Polynomial': g1, g2, g3, g4 = los_angeles(matrix([0.0, 0.0, 0.0, 0.15]), 'Polynomial')
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


def plot_results(delaytype):
    """Plot results from two_step_test()
    """
    
    if delaytype == 'Polynomial':
        # results for ind_obs=graph.indlinks.keys()
        tolls_collected1 = [[164.76921455176952, 690.919364787199, 1003.9715271443032], [164.76920846294715, 690.9196375266908, 1003.971503332828]]
        so_costs1 = [[1952.813985915695, 2991.422528084584, 4716.77478205995], [1952.813764068618, 2991.4225756924197, 4716.774942888717]]
        costs1 = [[1961.4088544768533, 3035.621746689017, 4751.282035585294], [1961.4088092322784, 3035.619283467948, 4751.282202058287]]
        ue_costs1 = [[1993.9477013089333, 3087.4591784440836, 4889.555908205202], [1993.9532544651581, 3087.4592286945117, 4889.562102247506]]
        
        # results for ind_obs=[(36,37,1), (13,14,1), (17,8,1), (24,17,1), (28,22,1), (14,13,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1)]
        tolls_collected2 = [[165.05685117271955, 354.1634370314711, 1015.9141858070569], [165.0567022829941, 354.1917223859498, 1015.9119098321756]]
        so_costs2 = [[1952.4102665477528, 2990.2968764530697, 4714.092736177128], [1952.813764068618, 2991.4225756924197, 4716.774942888717]]
        costs2 = [[1961.0182151023969, 3019.353015234613, 4750.051821659878], [1961.4090940730164, 3020.4668352547183, 4752.675436210774]]
        ue_costs2 = [[1993.5324364880598, 3086.401670312004, 4886.656461616093], [1993.9532544651581, 3087.4592286945117, 4889.562102247506]]
        
        # results for ind_obs=[(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
        tolls_collected3 = [[122.37203827953705, 369.5962596851687, 993.1260548913893], [124.00144138458633, 375.65796385498777, 1018.6147484597942]]
        so_costs3 = [[1819.720988976306, 2777.3939885366685, 4677.595258102356], [1952.813764068618, 2991.4225756924197, 4716.774942888717]]
        costs3 = [[1834.6351015726038, 2802.3016403566235, 4714.102145227426], [1965.8049083549738, 3020.126090534662, 4750.071391283308]]
        ue_costs3 = [[1872.3492340077084, 2897.9834368843217, 4886.769726452165], [1993.9532544651581, 3087.4592286945117, 4889.562102247506]]
        
        # results for ind_obs=[(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
        tolls_collected4 = [[602.0313412901077, 596.7300601881633, 833.3678073460065], [609.094482012133, 542.8579545671112, 843.1817520139325]]
        so_costs4 = [[13252.728980891832, 20178.927282098706, 28549.614451166774], [1952.813764068618, 2991.4225756924197, 4716.774942888717]]
        costs4 = [[13254.342842900576, 20239.30287102026, 28549.614731411177], [1992.3969007442322, 3130.0522527298062, 4805.6728306302175]]
        ue_costs4 = [[13259.433117552415, 20185.90937054816, 28556.779443740656], [1993.9532544651581, 3087.4592286945117, 4889.562102247506]]
        
    
    if delaytype == 'Hyperbolic':
        # results for ind_obs=graph.indlinks.keys() 
        tolls_collected1 = [[161.15742640370462, 416.80264009965265, 948.124688583883], [161.4013529462645, 417.25218843773354, 950.1674351402539]]
        so_costs1 = [[2519.5329969919685, 3679.4257321044365, 5325.940489651684], [2544.920874812245, 3718.292867749051, 5383.190501656683]]
        costs1 = [[2525.7590340322095, 3679.4261220632625, 5356.8071818038225], [2551.0113714206195, 3718.3028736461647, 5412.310857035818]]
        ue_costs1 = [[2548.8303153671195, 3736.67165067958, 5444.499035875823], [2574.155371287541, 3775.672784740404, 5503.564654849874]]
            
        # results for ind_obs=[(36,37,1), (13,14,1), (17,8,1), (24,17,1), (28,22,1), (14,13,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1)]
        tolls_collected2 = [[355.56391097372205, 701.3544774585599, 806.6640717434449], [374.34816990991186, 736.3763770517197, 842.0768966826593]]
        so_costs2 = [[2174.1387842065237, 3306.974167934956, 4969.5560792612605], [2544.920874812245, 3718.292867749051, 5383.190501656683]]
        costs2 = [[2201.273656230071, 3338.7536674213156, 4988.732654672769], [2564.0153714222406, 3746.115605602468, 5398.046518880896]]
        ue_costs2 = [[2209.622435372754, 3380.7186533697286, 5104.755602389524], [2574.155371287541, 3775.672784740404, 5503.564654849874]]
        
        # results for ind_obs=[(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
        tolls_collected3 = [[339.37838703976513, 737.9887571945254, 930.0531221139565], [353.28118085588176, 763.3708450755438, 955.644104623869]]
        so_costs3 = [[2520.1983680076073, 3833.234908857907, 5645.320052639795], [2544.920874812245, 3718.292867749051, 5383.190501656683]]
        costs3 = [[2536.138622203524, 3862.614483834193, 5668.430527845267], [2562.897499637454, 3750.9322754235955, 5406.402328441233]]
        ue_costs3 = [[2556.0266644367803, 3903.9758481975805, 5741.545038028959], [2574.155371287541, 3775.672784740404, 5503.564654849874]]
            
        # results for ind_obs=[(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
        tolls_collected4 = [[443.85669649193784, 837.5428841560781, 1380.117295219297], [479.57423528502, 931.2309202734045, 1499.6972922042573]]
        so_costs4 = [[2158.995768390798, 3843.4202195503713, 7595.092411648515], [2544.920874812245, 3718.292867749051, 5383.190501656683]]
        costs4 = [[2192.900868595389, 3888.0172745668465, 7626.17585256234], [2575.5748442305207, 3755.844017718783, 5450.304046277645]]
        ue_costs4 = [[2231.1674179937672, 3999.772800567613, 7867.415688946321], [2574.155371287541, 3775.672784740404, 5503.564654849874]]
    
    # effects of toll pricing
    opacity = 0.4
    width = 0.2
    index = np.arange(4)
    plt.xticks(index, ('all', '10', '4 good', '4 bad') )
    colors = ['y', 'r', 'b']
    demands = ['94', '118', '142']
    for i in [0,1,2]:
        so = so_costs1[1][i]
        ue = ue_costs1[1][i]
        delta = ue-so
        plt.bar(index-width/2+(i-1)*width,
                        [100.0*(costs1[1][i]-so)/delta, 
                        100.0*(costs2[1][i]-so)/delta, 
                        100.0*(costs3[1][i]-so)/delta, 
                        100.0*(costs4[1][i]-so)/delta],
                        width, alpha=opacity, color=colors[i], label='demand='+demands[i])
    plt.xlabel('Sets of observed links')
    plt.ylabel('Relative loss (%)')
    plt.title('Effects of toll pricing with '+delaytype+' delay function')
    plt.legend(loc=0)
    plt.show()
    
    # error in predicted UE
    for i in [0,1,2]:
        ue = ue_costs1[1][i]
        plt.bar(index-width/2+(i-1)*width,
                [100.0*abs(ue_costs1[0][i]-ue)/ue,
                 100.0*abs(ue_costs2[0][i]-ue)/ue,
                 100.0*abs(ue_costs3[0][i]-ue)/ue,
                 100.0*abs(ue_costs4[0][i]-ue)/ue],
                width, alpha=opacity, color=colors[i], label='demand='+demands[i])
    plt.xlabel('Sets of observed links')
    plt.ylabel('Relative error (%)')
    plt.title('Error in predicted UE with '+delaytype+' delay function')
    plt.legend(loc=0)
    plt.show()
    
    # error in predicted SO
    for i in [0,1,2]:
        so = so_costs1[1][i]
        plt.bar(index-width/2+(i-1)*width,
                [100.0*abs(so_costs1[0][i]-so)/so,
                 100.0*abs(so_costs2[0][i]-so)/so,
                 100.0*abs(so_costs3[0][i]-so)/so,
                 100.0*abs(so_costs4[0][i]-so)/so],
                width, alpha=opacity, color=colors[i], label='demand='+demands[i])
    plt.xlabel('Sets of observed links')
    plt.ylabel('Relative error (%)')
    plt.title('Error in predicted SO with '+delaytype+' delay function')
    plt.legend(loc=0)
    plt.show()
    
    # error in predicted tolled ue
    for i in [0,1,2]:
        plt.bar(index-width/2+(i-1)*width,
                [100.0*abs(costs1[0][i]-costs1[1][i])/costs1[1][i],
                 100.0*abs(costs2[0][i]-costs2[1][i])/costs2[1][i],
                 100.0*abs(costs3[0][i]-costs3[1][i])/costs3[1][i],
                 100.0*abs(costs4[0][i]-costs4[1][i])/costs4[1][i]],
                width, alpha=opacity, color=colors[i], label='demand='+demands[i])
    plt.xlabel('Sets of observed links')
    plt.ylabel('Relative error (%)')
    plt.title('Error in predicted tooled UE with '+delaytype+' delay function')
    plt.legend(loc=0)
    plt.show()


def main():
    theta = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
    graph = los_angeles()[3]
    ind_obs = {}
    ind_obs[0] = graph.indlinks.keys()
    ind_obs[1] = [(36,37,1), (13,14,1), (17,8,1), (24,17,1), (28,22,1), (14,13,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    ind_obs[2] = [(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    ind_obs[3] = [(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
    
    test_toll_pricing(theta, True)
    #two_step_test(ind_obs[3], 'Polynomial')
    #plot_results('Hyperbolic')
    
    


if __name__ == '__main__':
    main()