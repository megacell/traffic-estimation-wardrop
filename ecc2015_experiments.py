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
import toll_pricing as tp


a, b = 3.5, 3.0
coef = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
graph = los_angeles(coef, 'Polynomial')[0]


def experiment_estimation(indlinks_obs, delaytype, degree):
    """Demonstrate multi-objective optimization for structural estimation
    1. Generate the graph of L.A.
    2. Generate synthetic data in UE
    3. Apply multi-objective solver
    
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
    w_multi = [0.001, .01, .1, .5, .9, .99, 0.999] # weight on the observation residual
    r_gap, r_obs, x_est = invopt.multi_objective_solver([g1,g2,g3,g4], [x1[obs],x2[obs],x3[obs],x4[obs]], obs, degree, w_multi)
    print r_gap
    print r_obs
    u = matrix([x1,x2,x3,x4])
    r_est = [np.linalg.norm(u-x, 1) / np.linalg.norm(u, 1) for x in x_est]
    print r_est
    
    
def experiment_toll_pricing(ws_so, ws_toll):
    """Demonstrate multi-objective optimization for toll pricing model
    1. Generate the graph of L.A.
    2. Run multi-objective solver
    3. Post-process results
    4. Print results
    
    Parameters
    ----------
    ws_so: list of weights for so objective
    ws_toll: list of weight for toll objective
    """
    graph = los_angeles(coef, 'Polynomial')[3]
    r_gap, toll_est, loss_est, toll_res, loss_res = tp.multi_objective_solver(graph, coef, ws_so, ws_toll)
    for i in range(len(ws_so)):
        for j in range(len(ws_toll)):
            if r_gap[i,j] < 0.0: r_gap[i,j]=0.0
            if toll_est[i,j] < 0.0: toll_est[i,j]=0.0
            if loss_est[i,j] < 0.0: loss_est[i,j]=0.0
            if toll_res[i,j] < 0.0: toll_res[i,j]=0.0
            if loss_res[i,j] < 0.0: loss_res[i,j]=0.0
    print r_gap
    print toll_est
    print loss_est
    print toll_res
    print loss_res
    
    
def draw_tolls():
    """Draw tolls"""
    
    
def results_est_poly(degree):
    """Display results when the true delay is polynomial"""
    if degree==6:
        r_gap = [0.3514945478659001, 0.1961611531023568, 0.1716444145948751, 0.17065344166251223, 0.17045516975380648, 0.17048860117939105, 0.17097031561001252]
        r_obs = [0.59265598048104329, 0.033896958661912446, 0.00037350125226116084, 5.2242920696447353e-06, 1.5083106412345355e-07, 9.5201385120293674e-08, 1.0040413042368986e-07]
        r_est = [0.13937361069687632, 0.029465753080638112, 0.00034691238632399542, 0.00092494014720187403, 0.0011425727298296376, 0.0011396367098416948, 0.00071672673008782472]
    """
    if degree==5:
        r_gap = [0.49332717517052505, 0.26663569286718664, 0.23548597091235388, 0.2344311309136928, 0.2341984841323786, 0.23417456989816432, 0.23417166680733056]
        r_obs = [0.59261387804666921, 0.033839942215393444, 0.00037334986647596263, 5.2076372971657296e-06, 6.4654609106342581e-08, 5.3487448435150483e-10, 5.1460506905051785e-12]
        r_est = [0.13931300064229757, 0.029490373590458453, 0.00032884752988508178, 0.00091766810963546874, 0.0011615324432562722, 0.0011692145553356995, 0.0011872866939422449]
    
    if degree==4:
        r_gap = [0.5729653192001556, 0.32472035921088266, 0.29213100144782905, 0.2909219379299841, 0.2906801067552578, 0.2906545176536313, 0.2906515127309753]
        r_obs = [0.59260816465319988, 0.033813447416166639, 0.00037335606992238899, 5.2076427693730981e-06, 6.4655000887394377e-08, 5.3488142104726901e-10, 5.1461278636759956e-12]
        r_est = [0.13931435727239652, 0.029511197645028221, 0.00035143283454960644, 0.0009290440468794417, 0.0011614750902444535, 0.0011690849050645354, 0.0011871728103361389]
    
    if degree==3:
        r_gap = [0.6624564908797662, 0.4137402254189099, 0.3851899005249685, 0.3835720430596299, 0.3833099652246711, 0.38328161223874074, 0.38327896911927667]
        r_obs = [0.59260620929781926, 0.035161894883601046, 0.00037561798879280579, 5.1339230369548893e-06, 6.381497135990285e-08, 5.2881863320343698e-10, 5.2015732631997571e-12]
        r_est = [0.13931307665090023, 0.05165360310354588, 0.034953071312759668, 0.033971179436959528, 0.03382500546598248, 0.033809992698478418, 0.033807593991076919]
    
    if degree==2:
        r_gap = [0.26886751773390005, 0.16023857541850223, 0.1328394267284071, 0.13026443302848656, 0.13018781026840226, 0.13017895379745173, 0.13018237915979142]
        r_obs = [0.59260640769178774, 0.039012956159923999, 0.00027926365886101287, 3.4571895937728091e-06, 4.2954520920081236e-08, 3.6161144562592957e-10, 3.6266073277042652e-12]
        r_est = [0.1393131702138552, 0.077043817781165147, 0.06656234084299141, 0.066798748966516899, 0.066810819958887616, 0.06675085204849307, 0.066775816452977949]
    
    if degree==1:
        r_gap = [0.3600603622525474, 0.07398822061038958, 0.0629608973453844, 0.062485537172215826, 0.06247458817535771, 0.062473238185037354, 0.06247319387729183]
        r_obs = [0.67573329964465945, 0.013194910881419174, 8.9279582021092668e-05, 1.1124430182504824e-06, 1.381455724650183e-08, 1.1653459570704158e-10, 1.2436293540089436e-12]
        r_est = [0.1636889960092226, 0.14262366445679853, 0.20975780359601365, 0.21485693561169467, 0.21495011494538274, 0.21495737632039419, 0.21497898895670292]
    """
    
    
def results_est_hyper(degree):
    """Display results when the true delay is hyperbolic"""
    if degree==6:
        r_gap = [0.28131048379615914, 0.19360668630015063, 0.18294324324316735, 0.18135181038793996, 0.18130181867342157, 0.18129806279354277, 0.18133134602833872]
        r_obs = [0.41862638077775782, 0.013212508560595275, 0.00012708736874764273, 1.5688928281663424e-06, 8.1771883584478209e-08, 7.2836247584159581e-08, 9.3129623530239088e-08]
        r_est = [0.078284347401316984, 0.035892084227933531, 0.039703818272014593, 0.03940228808866747, 0.039395305610507671, 0.039397204940010318, 0.039397585861242083]
    """
    if degree==5:
        r_gap = [0.3928718016343673, 0.2683244905209487, 0.25315353254331346, 0.2511894374466198, 0.25112514576450307, 0.251117783870057, 0.25111786068339664]
        r_obs = [0.41848497751365027, 0.01320176986220244, 0.00012728784653201475, 1.4211775814973044e-06, 1.7626391929271559e-08, 1.4663541626984021e-10, 1.5088483641315054e-12]
        r_est = [0.078272311098310232, 0.035829949570359503, 0.039703375703081453, 0.039396016236798996, 0.039400213337798973, 0.039399739267838317, 0.039399832135254866]
    
    if degree==4:
        r_gap = [0.46285064780204327, 0.3284698437485553, 0.3112545333809956, 0.31069301337869126, 0.31062530793144205, 0.31061732719921975, 0.31061792962077855]
        r_obs = [0.41856629982569671, 0.013226346003668831, 0.00011555016032740922, 1.4498834449085235e-06, 1.7984788008509065e-08, 1.4976704560313976e-10, 1.5453790316245717e-12]
        r_est = [0.078325328989525975, 0.0357899476942923, 0.039553316704928135, 0.039438796834406896, 0.039436142461418931, 0.039435302193824434, 0.039437187467939727]
    
    if degree==3:
        r_gap = [0.5429875294781642, 0.4004771203995994, 0.38068248512879166, 0.3801095992202417, 0.38004582424294997, 0.3800384579273996, 0.3800392970034846]
        r_obs = [0.41852101043057971, 0.013231095733459633, 0.00011555090516328244, 1.4498837525376358e-06, 1.7984783121939954e-08, 1.4976720293600978e-10, 1.5453802841386389e-12]
        r_est = [0.078284224126328072, 0.035453003642433777, 0.039557824231840732, 0.039438716404762264, 0.039430813205593232, 0.03943722128227984, 0.039435130892829337]
    
    if degree==2:
        r_gap = [0.21466085418018863, 0.16731598937590406, 0.15707411700239557, 0.15669988150029024, 0.15666110493440719, 0.15665670575493423, 0.156656137719509]
        r_obs = [0.41843237946527051, 0.012968122748051045, 0.000120941813052212, 1.5141706652535714e-06, 1.8724200627051487e-08, 1.5764244004178764e-10, 1.5835812960406917e-12]
        r_est = [0.078570187077445133, 0.038964677742305602, 0.027914181466683372, 0.027527128025774544, 0.027480315329471566, 0.027475713118449915, 0.027475038815409861]

    
    if degree==1:
        r_gap = [0.2473728686809255, 0.14096147360607036, 0.13971982540757377, 0.13961205214936165, 0.1396001333465915, 0.13959882250087846, 0.13959992291384615]
        r_obs = [0.37468249556523475, 0.0068366420085976716, 8.0713371389038932e-05, 1.0379113008469728e-06, 1.2924891933922978e-08, 1.0844837983792924e-10, 1.173409614854846e-12]
        r_est = [0.090273325183187661, 0.060073922279314483, 0.059900322965394943, 0.059864074673187509, 0.059884762753935047, 0.059862465250625031, 0.059862577981822276]
    """


def results_toll():
    """Display results of experiment_toll_pricing
    
    optimal weights: w_so=1e2, w_toll=1e-2
    """
    r_gap = matrix([[ 0.00e+00,  4.32e-09,  3.01e-02,  9.82e-01,  9.65e-01],
    [ 0.00e+00,  4.34e-09,  3.01e-02,  3.20e-02,  9.65e-01],
    [ 0.00e+00,  7.02e-09,  2.58e-02,  3.20e-02,  9.65e-01],
    [ 3.17e-09,  4.34e-09,  3.08e-02,  3.21e-02,  9.65e-01],
    [ 6.43e-08,  4.44e-09,  3.08e-02,  3.21e-02,  9.65e-01]])
    
    toll_est = matrix([[ 7.89e+02,  7.90e+02,  5.19e+00,  0.00e+00,  0.00e+00],
    [ 9.22e+02,  7.89e+02,  5.19e+00,  0.00e+00,  0.00e+00],
    [ 7.31e+02,  7.31e+02,  3.45e+01,  0.00e+00,  0.00e+00],
    [ 6.79e+02,  6.79e+02,  4.46e+00,  0.00e+00,  0.00e+00],
    [ 6.78e+02,  6.78e+02,  4.45e+00,  0.00e+00,  0.00e+00]])
    
    loss_est = matrix([[ 1.46e-01,  1.49e-01,  1.12e-02,  2.92e-05,  7.45e-09],
    [ 2.04e-01,  1.46e-01,  1.09e-02,  1.55e-07,  7.28e-09],
    [ 3.39e-02,  3.46e-02,  3.61e-02,  3.36e-07,  1.49e-09],
    [ 1.13e-05,  1.15e-05,  0.00e+00,  0.00e+00,  5.31e-12],
    [ 0.00e+00,  1.18e-09,  0.00e+00,  0.00e+00,  4.81e-12]])
    
    toll_res = matrix([[ 7.89e+02,  7.90e+02,  3.41e+00,  9.68e-01,  0.00e+00],
    [ 9.22e+02,  7.89e+02,  3.40e+00,  0.00e+00,  0.00e+00],
    [ 7.31e+02,  7.31e+02,  3.57e+01,  0.00e+00,  0.00e+00],
    [ 6.79e+02,  6.79e+02,  2.61e+00,  0.00e+00,  0.00e+00],
    [ 6.78e+02,  6.78e+02,  2.60e+00,  0.00e+00,  0.00e+00]])
    
    loss_res = matrix([[ 1.49e-01,  1.51e-01,  9.95e-01,  1.00e+00,  1.00e+00],
    [ 2.04e-01,  1.48e-01,  9.95e-01,  1.00e+00,  1.00e+00],
    [ 3.43e-02,  3.57e-02,  8.11e-01,  1.00e+00,  1.00e+00],
    [ 1.19e-03,  3.29e-04,  9.91e-01,  1.00e+00,  1.00e+00],
    [ 1.20e-03,  1.29e-03,  9.91e-01,  1.00e+00,  1.00e+00]])
    
    log_r_gap = matrix(0.0, (5,5))
    for i in range(5):
        for j in range(5):
            log_r_gap[i,j]= max(-10.0, np.log10(r_gap[i,j]))
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.imshow(np.flipud(log_r_gap[:4,:4].T), interpolation='nearest', cmap=plt.cm.YlOrRd)
    plt.xlabel('Weight on toll minimization')
    plt.ylabel('Weight on SO objective')
    plt.colorbar()
    plt.title('Gap residual')
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.imshow(np.flipud(toll_est[:4,:4].T), interpolation='nearest', cmap=plt.cm.YlOrRd)
    plt.xlabel('Weight on toll minimization')
    plt.ylabel('Weight on SO objective')
    plt.colorbar()
    plt.title('Estimated toll')
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.imshow(np.flipud(loss_est[:4,:4].T), interpolation='nearest', cmap=plt.cm.YlOrRd)
    plt.xlabel('Weight on toll minimization')
    plt.ylabel('Weight on SO objective')
    plt.colorbar()
    plt.title('Estimated loss')
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.imshow(np.flipud(toll_res[:4,:4].T), interpolation='nearest', cmap=plt.cm.YlOrRd)
    plt.xlabel('Weight on toll minimization')
    plt.ylabel('Weight on SO objective')
    plt.colorbar()
    plt.title('Real toll')
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.imshow(np.flipud(loss_res[:4,:4].T), interpolation='nearest', cmap=plt.cm.YlOrRd)
    plt.xlabel('Weight on toll minimization')
    plt.ylabel('Weight on SO objective')
    plt.colorbar()
    plt.title('Real loss')
    plt.show()
    
    print r_gap
    print toll_est
    print loss_est
    print toll_res
    print loss_res

    
def main():
    #type = 'Hyperbolic'
    type = 'Polynomial'
    degree = 6
    ind_obs = [(36,37,1), (13,14,1), (17,8,1), (24,17,1), (28,22,1), (14,13,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    #ind_obs = [(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    #ind_obs = [(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
    #experiment_estimation(ind_obs, type, degree)
    #display_results_polynomial(degree)
    #display_results_hyperbolic(degree)
    experiment_toll_pricing([1e-4, 1e-2, 1e0, 1e2], [1e-4, 1e-2, 1e0, 1e2])
    #results_toll()


if __name__ == '__main__':
    main()