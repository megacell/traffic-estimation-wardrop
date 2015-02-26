'''
Created on Aug 10, 2014

@author: jeromethai
'''

"""The goal is to compare different path solvers
Generate synthetic data:
1. compute the UE link flow using node-link flow formulation
2. compute the 3-shortest paths for each OD pair with weights the delays in UE (optimal)
3. compute the UE path flow using link-path formulation with the paths enumerated above
4. check that the resulting link flow vectors are the same (yes!)
5. generate a trajectory of waypoints for each used path by taking the closest waypoints to each link along the path
6. gathering routes with same waypoints sequence and sum up the flow for these routes

We suppose that we have:
1. Topology of the network, #lanes, #ff_delays, OD pairs, OD flows (OSM, google maps, census data, tomogravity)
2. Network is in steady-state equilibrium (UE) (common assumption in traffic assignment)
3. We know the latency functions (inverse equilibrium)
4. We know exactly the used routes (google maps, K-shortest paths on UE delays)
5. We have waypoint trajectories and flows and the set of routes associated to each
6. We have partial link flows from static sensors

How can we estimate the path flows? 4 experiments:
1. min ||P*x-l|| s.t. U*x=r, x>=0 with U OD-path incidence (most under-determined)
2. Compute link-flow UE, then compute feasible path flows (very under-determined)
3. Compute UE path flows directly (very under-determined)
4. min ||P*x-l|| s.t. U*x=r, x>=0 with U WP-path incidence (weakly to non under-determined)

Remarks:
1. UE doesn't enable to determine path flows
2. Turning ratios not accessible
3. Hence path flows is not really accessible (few route estimation in the literature)
4. That's why 'route flow estimation' is not used in the literature
"""

#import ipdb
import matplotlib.pyplot as plt
import numpy as np
from cvxopt import matrix
import random

import Waypoints as WP
import wp_generator as wp
import ue_solver as ue
import path_solver as path
import rank_nullspace as rn
from generate_graph import los_angeles
from generate_paths import find_UESOpaths
from path_solver import linkpath_incidence

data = []
#data.append((30, 60, 0.2, [((3.5, 0.5, 6.5, 3.0), 30)], (12,6), 2.0))
#data.append((20, 40, 0.2, [((3.5, 0.5, 6.5, 3.0), 20)], (10,5), 2.0))
data.append((10, 20, 0.2, [((3.5, 0.5, 6.5, 3.0), 10)], (8,4), 2.0))
data.append((5, 10, 0.2, [((3.5, 0.5, 6.5, 3.0), 5)], (4,2), 2.0))
data.append((3, 5, 0.2, [((3.5, 0.5, 6.5, 3.0), 2)], (4,2), 2.0))
# data[4] = (1, 3, 0.2, [((3.5, 0.5, 6.5, 3.0), 1)], (2,2), 2.0)
theta = matrix([0.0, 0.0, 0.0, 0.15])


def synthetic_data(data=None, SO=False, demand=3, N=10, path=None):
    """Generate synthetic data for the experiments
    
    Parameters:
    -----------
    data: (N0, N1, scale, regions, res, margin)
        N0: number of background samples
        N1: number of samples on links
        regions: list of regions, regions[k] = (geometry, N_region)
        res: (n1, n2) s.t. the width is divided into n1 cells and the height into n2 cells
        margin: margin around each cell
    SO: if True, computes SO
    demand: OD demand
    N: number of sets of observed links
    e.g. if N=10, we run 10 experiments with observations from the 10, 20, ..., 100% most occupied links
    
    Return value:
    ------------
    g: Graph object with paths
    f_true: true path flow
    l_true: UE linkflow
    path_wps: dictionary {path_id: wp_ids}, with wp_ids list of waypoints
    wp_trajs: waypoint trajectories [(wp_traj, path_list, flow)]
    obs: dictionary {i: [link_ids]} where obs[i] are the indices of the i*n/N links with the most flow
    """
    random = True
    g, f_true, path_wps, wp_trajs, l_true = wp.compute_wp_flow(SO, demand, random, data, path=path)
    n = g.numlinks
    g.visualize()
    # P = linkpath_incidence(g)
    # l = P*f_true  # for numerical error purposes
    obs, sorted = {}, np.array(l_true).argsort(axis=0)
    for i in range(N-1):
        tmp = np.sort(sorted[-(i+1)*n/N:], axis=0)[:,0]
        obs[i] = [int(k) for k in tmp]
    obs[N-1] = range(n)
    return g, f_true, l_true, path_wps, wp_trajs, obs


def ratio_wptrajs_usedpaths(trials=10, demand=3):
    ratiosUE, ratiosSO = [0.0]*len(data), [0.0]*len(data)
    for k in range(trials):
        for i in range(5):
            g, f_true, path_wps, wp_trajs = wp.compute_wp_flow(False, demand, True, data[i])
            ratiosUE[i] += float(100*len(wp_trajs)) / len(path_wps)
            g, f_true, path_wps, wp_trajs = wp.compute_wp_flow(True, demand, True, data[i])
            ratiosSO[i] += float(100*len(wp_trajs)) / len(path_wps)
    for i in range(5):
        ratiosUE[i] /= trials
        ratiosSO[i] /= trials
    print ratiosUE
    print ratiosSO

def experiment(data=None, SO=False, trials=5, demand=3, N=10, withODs=False, data_id=None):
    """Run set of experiments
    Steps:
    1. generate synthetic data with synthetic_data()
        g: Graph object with paths
        f_true: true path flow
        l: UE linkflow
        path_wps: dictionary {path_id: wp_ids}, with wp_ids list of waypoints
        wp_trajs: waypoint trajectories [(wp_traj, path_list, flow)]
        obs: dictionary {i: [link_ids]} where obs[i] are the indices of the i*n/N links with the most flow
    2. Solve the linear regression problem min ||P*f-l|| s.t. U*f=r, x>=0
        a. solve it with {i*n/N, i=1,...,N} observed links and OD information
        b. solve it with {i*n/N, i=1,...,N} observed links and cell-path flow (+ OD flows if withODs)
    3. Compute the relative link flow errors l_errors and the relative path flow errors
    4. Repeat this 'trials' (=10) times and compute the average error and standard deviation
    
    Parameters:
    ----------
    data: (N0, N1, scale, regions, res, margin) inputs of generate_wp
    SO: if True, computes SO
    trials: number of trials and take the average
    demand: OD demand
    
    N: number of sets of observed links
    e.g. if N=10, we run 10 experiments with observations from the 10, 20, ..., 100% most occupied links
    
    plot: if True, plot results
    withODs: include ODs in the constraints
    data_id: id when figure is saved
    """
    numexp = 2*N
    l_errors, f_errors = [[] for i in range(numexp)], [[] for i in range(numexp)]
    mean_ratio = 0
    ddl_ODs, ddl_cellpaths = [[] for i in range(N)], [[] for i in range(N)]
    k = 0
    while k < trials:
        print 'trial', k
        g, f_true, l, path_wps, wp_trajs, obs = synthetic_data(data, SO, demand, N)
        mean_ratio += float(len(wp_trajs)) / len(path_wps)
        norm_l, norm_f = np.linalg.norm(l, 1), np.linalg.norm(f_true, 1)
        err_f = lambda x: np.linalg.norm(f_true-x, 1) / norm_f
        err_l = lambda x: np.linalg.norm(l-x, 1) / norm_l
        P = path.linkpath_incidence(g)
        ls, fs = [], []
        #print 'Compute min ||P*f-l|| s.t. U*f=r, x>=0 with U=OD-path incidence matrix'
        # A trial must complete successfully for all 2*N tests for it to count
        failed = False
        for i in range(N):
            try:
                f, rank, dim = path.feasible_pathflows(g, l[obs[i]], obs[i], with_ODs=True, x_true=f_true)
                ddl_ODs[i].append(dim - rank)
            except (ValueError, UnboundLocalError) as e:
                print e
                # 'Probably your QP is either non-positively defined (for cvxopt_qp you should have xHx > 0 for any x != 0) or very ill-conditioned.'           # __str__ allows args to be printed directly
                failed = True
                break
            fs.append(f)
            ls.append(P*f)
        #print 'Compute min ||P*f-l|| s.t. U*f=r, x>=0 with U=waypoint-path incidence matrix'
        if failed:
            continue
        for i in range(N):
            try:
                f, rank, dim = path.feasible_pathflows(g, l[obs[i]], obs[i], 
                    with_ODs=withODs, with_cell_paths=True, x_true=f_true, wp_trajs=wp_trajs)
                ddl_cellpaths[i].append(dim - rank)
                # Throw out trials for which any domain error is caught
                #for j in range(U.size[0]):
                #    path.feasible_pathflows(g, l[obs[i]], obs[i], eq_constraints=(U[:j,:],r[:j]))
            except (ValueError, UnboundLocalError) as e:
                print e
                # 'Probably your QP is either non-positively defined (for cvxopt_qp you should have xHx > 0 for any x != 0) or very ill-conditioned.'           # __str__ allows args to be printed directly
                failed = True
                break
            fs.append(f)
            ls.append(P*f)

        if failed:
            continue
        k += 1
        l_error = [err_l(ls[i]) for i in range(numexp)]
        f_error = [err_f(fs[i]) for i in range(numexp)]
        """
        if withODs:
            print f_error[:N], f_error[N:]
            assert np.all([bo<=od for od,bo in zip(f_error[:N],f_error[N:])]), \
                'Configurations with both OD+cellpath should be performing ' + \
                'at least as well as with just OD'
            assert np.all([bo<=od for od,bo in zip(l_error[:N],l_error[N:])]), \
                'Configurations with both OD+cellpath should be performing ' + \
                'at least as well as with just OD'
        """
        for i in range(numexp):
            l_errors[i].append(l_error[i])
            f_errors[i].append(f_error[i])
    mean_l_errors = [np.mean(l_errors[i]) for i in range(numexp)]
    mean_f_errors = [np.mean(f_errors[i]) for i in range(numexp)]
    std_l_errors = [np.std(l_errors[i]) for i in range(numexp)]
    std_f_errors = [np.std(f_errors[i]) for i in range(numexp)]
    mean_ddl_ODs = [np.mean(ddl_ODs[i]) for i in range(numexp)]
    mean_ddl_cellpaths = [np.mean(ddl_cellpaths[i]) for i in range(numexp)]
    return mean_l_errors, mean_f_errors, std_l_errors, std_f_errors,\
        mean_ratio/trials, mean_ddl_ODs, mean_ddl_cellpaths


def run_experiments(SO=False, trials=5):
    """
    Run experiment(data=data[i]), i=0,...,4 with decreasing densities of waypoints
    
    Output:
    -------
    A: mean linkflow errors over 10 trials
       A[i,j-1] for data[i], j*n/N observed links, j = 1,...,10 with ODs alone
       A[i,10+j-1] for data[i], j*n/N observed links, j = 1,...,10 with cellpath flows
    B: mean pathflow errors over 10 trials
       B[i,j-1] for data[i], j*n/N observed links, j = 1,...,10 with ODs alone
       B[i,10+j-1] for data[i], j*n/N observed links, j = 1,...,10 with cellpath flows
    C: std linkflow errors over 10 trials
       C[i,j-1] for data[i], j*n/N observed links, j = 1,...,10 with ODs alone
       C[i,10+j-1]for data[i], j*n/N observed links, j = 1,...,10 with cellpath flows
    D: std pathflow errors over 10 trials
       D[i,j-1] for data[i], j*n/N observed links, j = 1,...,10 with ODs alone
       D[i,10+j-1] for data[i], j*n/N observed links, j = 1,...,10 with cellpath flows
    n is the number of links
    N number of sets of observed links
    e.g. if N=10, we run 10 experiments with observations from the 10, 20, ..., 100% most occupied links
    """
    if SO: ue_or_so = 'SO'
    else: ue_or_so = 'UE'
   
    with open('ISTTT_results/ISTTT_results_' + ue_or_so + '.txt', 'w') as out:
        A, B, C, D, E, F, G = [], [], [], [], [], [], []
        for i in range(len(data)):
            print 'Experiment with data[{}] without ODs'.format(i)
            a, b, c, d, e, f, g = experiment(data[i], SO=SO, withODs=False, data_id=i, trials=trials)
            A.append(a); B.append(b); C.append(c); D.append(d), E.append(e), F.append(f), G.append(g)
            print 'Experiment with data[{}] with ODs'.format(i)
            a, b, c, d, e, f, g = experiment(data[i], SO=SO, withODs=True, data_id=i, trials=trials)
            A.append(a); B.append(b); C.append(c); D.append(d), E.append(e), F.append(f), G.append(g)
        # Write data to output
        for a in A: out.write('%s\n' % a)  # mean
        for b in B: out.write('%s\n' % b)  # mean
        for c in C: print c  # standard deviation
        for d in D: print d  # standard deviation
        
    with open('ISTTT_results/ISTTT_ddl_ODs_' + ue_or_so + '.txt', 'w') as out:
        for f in F: out.write('%s\n' % f)
    with open('ISTTT_results/ISTTT_ddl_cellpaths_' + ue_or_so + '.txt', 'w') as out:
        for g in G: out.write('%s\n' % g)
    with open('ISTTT_results/ISTTT_ratio_' + ue_or_so + '.txt', 'w') as out:
        out.write('%s' % E)


def QP_rank(data=None, SO=False, trials=10, demand=3, N=10, withODs=False, wp_model=True):
    """Compute the rank of the QP used to solve the problem with or without waypoint data
    
    Parameters:
    -----------
    data: (N0, N1, scale, regions, res, margin) inputs of generate_wp
    SO: if True, computes SO
    trials: number of trials and take the average
    demand: OD demand
    N: number of bins
    
    Return value:
    -------------
    rank: rank of the QP
    dim: dimension of the problem
    """
    ranks = [0]*N
    if wp_model:
        num_used_paths = 0
        for k in range(trials):
            g, f_true, l, path_wps, wp_trajs, obs = synthetic_data(data, SO, demand, N)
            P = path.linkpath_incidence(g)
            if k==0: dim = P.size[1]
            U,r = WP.simplex(g, wp_trajs, withODs)
            for i in range(N): ranks[i] += rn.rank(matrix([P[obs[i],:], U]))
            num_used_paths += len(path_wps)
        for i in range(N): ranks[i] = float(ranks[i])/trials
        num_used_paths = float(num_used_paths)/trials
        return ranks, dim, num_used_paths
    else:
        g = los_angeles(theta, 'Polynomial')[demand]
        paths = find_UESOpaths(SO)
        for p in paths: g.add_path_from_nodes(p)
        l = ue.solver(g, update=True, SO=SO)
        n = g.numlinks
        P = path.linkpath_incidence(g)
        U,r = path.simplex(g)
        dim = P.size[1]
        obs, sorted = {}, np.array(l).argsort(axis=0)
        for i in range(N-1):
            tmp = np.sort(sorted[-(i+1)*n/N:], axis=0)[:,0]
            obs[i] = [int(k) for k in tmp]
        obs[N-1] = range(n)
        for i in range(N): ranks[i] = rn.rank(matrix([P[obs[i],:], U]))
        return ranks, dim
        
        
def run_QP_ranks(wp_model=True):
    """Compute the ranks, dimensions, and number of used paths
    for the QP problem with waypoints 
    """
    SO = [False, False, True, True]
    ODs = [False, True, False, True]
    ranks, dims, num_used_paths = [], [], []
    for i in range(5):
        for k in range(4):
            if wp_model: rank, dim, n = QP_rank(data[i], SO[k], withODs=ODs[k], wp_model=wp_model)
            else: rank, dim = QP_rank(data[i], SO[k], withODs=ODs[k], wp_model=wp_model)
            ranks.append(rank)
            dims.append(dim)
            if wp_model: num_used_paths.append(n)
    print ranks
    print dims
    print num_used_paths


def display_ranks():
    num_links, num_ods = 122, 42
    num_obs = [13, 25, 37, 49, 61, 74, 86, 98, 110, 122]
    index, color, num_wps = [10*i for i in range(1,11)], ['m', 'c', 'b', 'k', 'g'], [80, 40, 20, 10, 5]
    state = ['UE', 'UE', 'SO', 'SO']
    ODs = ['without', 'with', 'without', 'with']
    ranks = [[101.8, 113.8, 124.8, 131.8, 141.8, 148.8, 151.8, 153.8, 154.8, 154.8], [138.0, 148.0, 157.0, 162.0, 165.0, 168.0, 170.0, 170.0, 170.0, 170.0], [163.3, 174.3, 182.3, 188.3, 195.3, 204.3, 205.3, 207.3, 208.3, 208.3], [184.5, 194.5, 201.7, 205.7, 209.7, 212.7, 212.7, 212.7, 212.7, 212.7], [95.2, 107.2, 118.2, 125.2, 135.4, 142.4, 145.4, 147.4, 148.4, 148.4], [128.5, 138.5, 147.5, 152.5, 155.5, 158.5, 160.5, 160.5, 160.5, 160.5], [141.0, 152.1, 160.6, 166.7, 173.7, 182.8, 183.8, 185.8, 186.8, 186.8], [171.8, 182.0, 189.2, 193.2, 197.2, 200.2, 200.2, 200.2, 200.2, 200.2], [77.7, 89.7, 100.7, 107.7, 118.1, 125.1, 128.1, 130.1, 131.1, 131.1], [113.2, 123.2, 132.5, 137.5, 140.5, 143.5, 145.5, 145.5, 145.5, 145.5], [99.4, 111.1, 119.8, 126.4, 133.4, 142.8, 143.8, 145.8, 146.8, 146.8], [140.7, 151.0, 158.7, 163.1, 167.1, 170.1, 170.1, 170.1, 170.1, 170.1], [56.2, 68.2, 79.7, 86.7, 97.3, 104.3, 107.3, 109.3, 110.3, 110.3], [91.8, 101.8, 111.4, 116.4, 119.4, 122.4, 124.4, 124.4, 124.4, 124.4], [64.7, 76.6, 86.5, 93.9, 100.9, 110.6, 111.6, 113.6, 114.6, 114.6], [90.9, 101.7, 110.7, 115.2, 119.2, 122.4, 122.4, 122.4, 122.4, 122.4], [28.0, 40.0, 51.9, 58.9, 69.8, 76.8, 79.8, 81.8, 82.8, 82.8], [70.6, 80.6, 91.0, 96.0, 99.0, 102.0, 104.0, 104.0, 104.0, 104.0], [36.2, 48.2, 59.1, 67.1, 74.1, 84.7, 85.7, 87.7, 88.7, 88.7], [74.3, 85.1, 94.5, 99.5, 103.5, 106.8, 106.8, 106.8, 106.8, 106.8]]
    ranks2 = [[52, 62, 73, 78, 81, 84, 86, 86, 86, 86], [52, 62, 73, 78, 81, 84, 86, 86, 86, 86], [50, 61, 71, 77, 81, 85, 85, 85, 85, 85], [50, 61, 71, 77, 81, 85, 85, 85, 85, 85], [52, 62, 73, 78, 81, 84, 86, 86, 86, 86], [52, 62, 73, 78, 81, 84, 86, 86, 86, 86], [50, 61, 71, 77, 81, 85, 85, 85, 85, 85], [50, 61, 71, 77, 81, 85, 85, 85, 85, 85], [52, 62, 73, 78, 81, 84, 86, 86, 86, 86], [52, 62, 73, 78, 81, 84, 86, 86, 86, 86], [50, 61, 71, 77, 81, 85, 85, 85, 85, 85], [50, 61, 71, 77, 81, 85, 85, 85, 85, 85], [52, 62, 73, 78, 81, 84, 86, 86, 86, 86], [52, 62, 73, 78, 81, 84, 86, 86, 86, 86], [50, 61, 71, 77, 81, 85, 85, 85, 85, 85], [50, 61, 71, 77, 81, 85, 85, 85, 85, 85], [52, 62, 73, 78, 81, 84, 86, 86, 86, 86], [52, 62, 73, 78, 81, 84, 86, 86, 86, 86], [50, 61, 71, 77, 81, 85, 85, 85, 85, 85], [50, 61, 71, 77, 81, 85, 85, 85, 85, 85]]
    dims = [275, 275, 300, 300, 275, 275, 300, 300, 275, 275, 300, 300, 275, 275, 300, 300, 275, 275, 300, 300]
    num_used = [91.0, 91.0, 153.0, 153.0, 91.0, 91.0, 153.0, 153.0, 91.0, 91.0, 153.0, 153.0, 91.0, 91.0, 153.0, 153.0, 91.0, 91.0, 153.0, 153.0]
    for k in range(4):
        dim, used = dims[k], num_used[k]
        plt.plot(index, [0]*10, '--k')
        plt.plot(index, [round(used-j, 1) for j in ranks2[k]], '--r', label='With ODs')
        for i in range(5):
            data = [round(used-j, 1) for j in ranks[4*i+k]]
            plt.plot(index, data, '-o'+color[i], label='With {} cells'.format(num_wps[i]))
        plt.title('Network in '+state[k]+' '+ODs[k]+' ODs, dim: {}, used paths: {}'.format(dims[4*i+k], num_used[4*i+k]))
        plt.xlabel('Percentage of links observed (%)')
        plt.ylabel('Degrees of freedom')
        plt.legend(loc=0)
        plt.show()
    #print len(ranks)


def display_ratios():
    ratios1 = [96.59340659340658, 92.52747252747251, 75.27472527472527, 48.35164835164835, 21.318681318681318]
    ratios2 = [92.6797385620915, 83.33333333333334, 59.28104575163398, 38.104575163398685, 14.11764705882353]
    labels = ['5', '10', '20', '40', '80']
    index = range(5)
    plt.plot(index, ratios1[::-1], '-o', label='UE')
    plt.plot(index, ratios2[::-1], '-o', label='SO')
    plt.title('Number of Waypoints over number of used paths')
    plt.xlabel('Number of cells')
    plt.ylabel('Percentage')
    plt.xticks(index, labels)
    plt.legend(loc=0)
    plt.show()
    

def display_results():
    # FORMAT of output file
    # The first 'N=10' lines are link flow error
    # The second 'N=10' lines are route flow error
    # With each set of N lines, we have 2*len(data) lines, alternating between
    # results without ODs and with ODs
    # The first 10 entries are with OD flows only
    # The last 10 entries are with an equality constraint too
    D = len(data)
    index = [10*i for i in range(1,11)]
    color = ['m', 'c', 'b', 'k', 'g']
    num_wps = [d[0] + d[1] + d[3][0][1] for d in data]

    # collect results for UE-type behavior
    for mode in ['UE', 'SO']:
        
        results = open('ISTTT_results/ISTTT_results_' + mode + '.txt', 'r').readlines()[2*D:]
        results  = [[float(e) for e in r[1:-2].split(', ')] for r in results]
        est_wp_woODs = [results[2*i][10:] for i in range(D)] #
        est_wp_wiODs = [results[2*i+1][10:] for i in range(D)]
        est_lf = [results[2*i][:10] for i in range(D)]
        
        for est_wp, str in [(est_wp_woODs, 'cellpath'), (est_wp_wiODs, 'OD+cellpath')]:
            plt.plot(index, est_lf[0], '-or', label='With OD flows')
            for i in range(D):
                plt.plot(index, est_wp[i], '-o'+color[i], label='With {} cells'.format(num_wps[i]))
            plt.title('Path flow errors for network in ' + mode + ': OD vs ' + str)
            plt.xlabel('Percentage of links observed (%)')
            plt.ylabel('Relative error')
            plt.yscale('log')
            plt.legend(loc=0)
            plt.show()
    
    # collect ranks for UE-type behavior
    for mode in ['UE', 'SO']:
        ddls = open('ISTTT_results/ISTTT_ddl_ODs_' + mode + '.txt', 'r').readlines()[0]
        ddl_ODs = [float(ddl) for ddl in ddls[1:-2].split(', ')][:10]
        ddls = open('ISTTT_results/ISTTT_ddl_cellpaths_' + mode + '.txt', 'r').readlines()
        ddls = [[float(r) for r in ddl[1:-2].split(', ')][:10] for ddl in ddls]
        ddl_cellpaths_woODs = [ddls[2*i] for i in range(D)]
        ddl_cellpaths_wiODs = [ddls[2*i+1] for i in range(D)]
        
        for ddl_wp, str in [(ddl_cellpaths_woODs, 'cellpath'), (ddl_cellpaths_wiODs, 'OD+cellpath')]:
            plt.plot(index, [round(r) for r in ddl_ODs], '--r', label='With ODs')
            for i in range(D):
                plt.plot(index, [round(r) for r in ddl_wp[i]], '-o'+color[i], label='With {} cells'.format(num_wps[i]))
            plt.title('Degree of freedom for network in ' + mode + ': OD vs ' + str)
            plt.xlabel('Percentage of links observed (%)')
            plt.ylabel('Degrees of freedom')
            plt.legend(loc=0)
            plt.show()



def main():
    myseed=29347293
    np.random.seed(myseed)
    random.seed(myseed)

    trials = 5
    #synthetic_data()
    #experiment()
    #ratio_wptrajs_usedpaths()

    # ISTTT experiments:
    SO = True
    run_experiments(SO=SO, trials=trials)
    display_results()

    #run_QP_ranks(False)
    #display_ranks()
    #display_ratios()

if __name__ == '__main__':
    main()