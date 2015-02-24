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

import ipdb
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


data = []
data.append((30, 60, 0.2, [((3.5, 0.5, 6.5, 3.0), 30)], (12,6), 2.0))
data.append((20, 40, 0.2, [((3.5, 0.5, 6.5, 3.0), 20)], (10,5), 2.0))
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
    l: UE linkflow
    path_wps: dictionary {path_id: wp_ids}, with wp_ids list of waypoints
    wp_trajs: waypoint trajectories [(wp_traj, path_list, flow)]
    obs: dictionary {i: [link_ids]} where obs[i] are the indices of the i*n/N links with the most flow
    """
    random = True
    
    g, f_true, path_wps, wp_trajs = wp.compute_wp_flow(SO, demand, random, data, path=path)
    l = ue.solver(g, update=True, SO=SO)
    n = g.numlinks
    g.visualize()
    obs, sorted = {}, np.array(l).argsort(axis=0)
    for i in range(N-1):
        tmp = np.sort(sorted[-(i+1)*n/N:], axis=0)[:,0]
        obs[i] = [int(k) for k in tmp]
    obs[N-1] = range(n)
    return g, f_true, l, path_wps, wp_trajs, obs


def ratio_wptrajs_usedpaths(trials=10, demand=3):
    ratiosUE, ratiosSO = [0.0]*5, [0.0]*5
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

def experiment(data=None, SO=False, trials=5, demand=3, N=10, plot=False,
               withODs=False, data_id=None):
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
    k = 0
    while k < trials:
        print 'trial', k
        g, f_true, l, path_wps, wp_trajs, obs = synthetic_data(data, SO, demand, N)
        norm_l, norm_f = np.linalg.norm(l, 1), np.linalg.norm(f_true, 1)
        P = path.linkpath_incidence(g)
        U,r = WP.simplex(g, wp_trajs, withODs)
        ls, fs = [], []
        #print 'Compute min ||P*f-l|| s.t. U*f=r, x>=0 with U=OD-path incidence matrix'

        # A trial must complete successfully for all 2*N tests for it to count
        failed = False
        for i in range(N):
            try:
                f = path.feasible_pathflows(g, l[obs[i]], obs[i])
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
                f = path.feasible_pathflows(g, l[obs[i]], obs[i], eq_constraints=(U,r))
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
        for i in range(numexp):
            l_errors[i].append(np.linalg.norm(l-ls[i], 1) / norm_l)
            f_errors[i].append(np.linalg.norm(f_true-fs[i], 1) / norm_f)
    mean_l_errors = [np.mean(l_errors[i]) for i in range(numexp)]
    mean_f_errors = [np.mean(f_errors[i]) for i in range(numexp)]
    std_l_errors = [np.std(l_errors[i]) for i in range(numexp)]
    std_f_errors = [np.std(f_errors[i]) for i in range(numexp)]
    if plot: plot_results(mean_l_errors, mean_f_errors, std_l_errors, std_f_errors, numexp, SO, data_id)
    return mean_l_errors, mean_f_errors, std_l_errors, std_f_errors


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
    if SO == True:
        output = 'ISTTT_results/ISTTT_results_SO.txt'
    else:
        output = 'ISTTT_results/ISTTT_results_UE.txt'

    with open(output, 'w') as out:
        A, B, C, D = [], [], [], []
        for i in range(len(data)):
            print 'Experiment with data[{}] without ODs'.format(i)
            a, b, c, d = experiment(data[i], SO=SO, plot=False, withODs=False,
                                    data_id=i, trials=trials)
            A.append(a)
            B.append(b)
            C.append(c)
            D.append(d)
            print a
            print b
            print 'Experiment with data[{}] with ODs'.format(i)
            a, b, c, d = experiment(data[i], SO=SO, plot=False, withODs=True,
                                    data_id=i, trials=trials)
            A.append(a)
            B.append(b)
            C.append(c)
            D.append(d)
            print a
            print b
        # Write data to output
        for a in A:
            out.write('%s\n' % a)  # mean
        for b in B:
            out.write('%s\n' % b)  # mean
        for c in C:
            print c  # standard deviation
        for d in D:
            print d  # standard deviation


def plot_results(mean_l_errors, mean_f_errors, std_l_errors, std_f_errors, numexp, SO, id):
    """Plot results of the experiment"""
    if SO: state = 'SO'
    else: state = 'UE'
    k = 1
    for log in [True, False]:
        opacity = 0.4
        bar_width = 0.35
        N = numexp/2
        index = np.arange(N)
        plt.bar(index, mean_l_errors[:N], bar_width, alpha=opacity,
                     color='b', log=log, label='With OD demands')
        plt.bar(index+bar_width, mean_l_errors[N:], bar_width, alpha=opacity,
                     color='r', log=log, label='With waypoints')
        plt.title('Link flow errors for network in ' + state)
        plt.xlabel('Percentage of links observed (%)')
        plt.ylabel('Relative error')
        plt.xticks(index + bar_width, [str(100*(i+1)/N) for i in range(N)])
        plt.legend()
        plt.savefig('../Desktop/'+state+str(id)+str(k)+'.eps')
        plt.clf()
        k += 1
        
        plt.bar(index, mean_f_errors[:N], bar_width, alpha=opacity,
                     color='b', log=log, label='With OD demands')
        plt.bar(index+bar_width, mean_f_errors[N:], bar_width, alpha=opacity,
                     color='r', log=log, label='With waypoints')
        plt.title('Path flow errors for network in ' + state)
        plt.xlabel('Percentage of links observed (%)')
        plt.ylabel('Relative error')
        plt.xticks(index + bar_width, [str(100*(i+1)/N) for i in range(N)])
        plt.legend()
        plt.savefig('../Desktop/'+state+str(id)+str(k)+'.eps')
        plt.clf()
        k += 1


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



def display_results():
    wpUE, wpSO, numexp = [], [], 20
    
    list =[[1.1918609582876509, 1.0334145777655983, 0.53709960179416394, 0.34911945621545332, 0.24056328327550097, 0.086763879114630613, 0.024370326143836753, 0.022596956275376539, 0.02199816114086146, 0.021714136223052467, 0.001955755504758361, 0.0020891070261361178, 0.0019447231253520742, 0.0019160088671532772, 0.0010530196135859495, 0.0010370457889297743, 0.0010364726009813269, 0.0010390184924027978, 0.0010292763356023965, 0.0010287920353146392],
    [1.1481852192519812, 0.82543698874225702, 0.60626605899741948, 0.51458825532598851, 0.14609039144388744, 0.059628076657257147, 0.059258267892067215, 0.059093668041230826, 0.058989156719217249, 0.058180965417447425, 0.0078290920303092544, 0.0077818814996883133, 0.0069394316898837797, 0.0063562527592264748, 0.0070376984508603259, 0.0069900032717506982, 0.0069400918947932124, 0.0069497239473369516, 0.0065802708848285551, 0.0065823533939581343],
    [1.1912184208224528, 1.0326877227245155, 0.53670114140765379, 0.34869509250715847, 0.23954460507041517, 0.08561965970980176, 0.023654938780723948, 0.02182616688832173, 0.021131551314520446, 0.020778264722955637, 0.020306716515900495, 0.016553109550921742, 0.013503795526295562, 0.012847670733480362, 0.0065871949359553482, 0.0066533001225822855, 0.0066541316820041526, 0.0066568858810975061, 0.0066469500212518549, 0.0066420709384414104],
    [1.1483195514835853, 0.82628259643899449, 0.60632854727769325, 0.51517119063945527, 0.14671651697296623, 0.060565246826714617, 0.060203272065878131, 0.059972940627361762, 0.059870338733122509, 0.059064169738999937, 0.071246783203145217, 0.033698866504999522, 0.014067823057872519, 0.013042670273942051, 0.014108778185643697, 0.013051244073416643, 0.013011477015841452, 0.012554025231540214, 0.012533887206507038, 0.01257144049673615],
    [1.1915284938135042, 1.0330156080282218, 0.53690257956029852, 0.34892079428199274, 0.24011765101756288, 0.086156361921374089, 0.024017725119869986, 0.022178421114286268, 0.021530742187971567, 0.021137467019568087, 0.10190680163450117, 0.062327122188556294, 0.049489965367358786, 0.048972645643363656, 0.039664695110966632, 0.029622387800151008, 0.02960538396891222, 0.029700839458719635, 0.029731523215007789, 0.029723522349242037],
    [1.1484164526238239, 0.82528243909259147, 0.60645334831889064, 0.51463205594554218, 0.14600692940596147, 0.0597258032110905, 0.059345882769775649, 0.059178436605596796, 0.059077285822027828, 0.05826475163529906, 0.24211800469718972, 0.10641490062912454, 0.03909307726417692, 0.035836633532026659, 0.034662472271925281, 0.029298929666913821, 0.029205361355424609, 0.027505201926456623, 0.0284758194449658, 0.028421510628689993],
    [1.1917494587964292, 1.0330734597734055, 0.53685980104918263, 0.34916949891171922, 0.24019857697931304, 0.086368479927934061, 0.024515134895970021, 0.022669170457775018, 0.022004327700394625, 0.021575714249672261, 0.27988210697484583, 0.20041972905478392, 0.12205937178557447, 0.11484401433111262, 0.10014559910781984, 0.087381574228058451, 0.087377604908392165, 0.087479928915239369, 0.087506848488304728, 0.087538150141245802],
    [1.1483087677030059, 0.82533182898227719, 0.60636268805978288, 0.51456244984372768, 0.14581497471693117, 0.059492114381458829, 0.059104490290057786, 0.05894465503090425, 0.058842342409228456, 0.058013646805397355, 0.37696347240777095, 0.17529681601757083, 0.083342390132848776, 0.076665254101201744, 0.071110085862349889, 0.056349843786007867, 0.056083886532416513, 0.054470125532995309, 0.054899201874500848, 0.054824629455407926],
    [1.1912350235200364, 1.0326792696685747, 0.53649916208294324, 0.34858752810027782, 0.23942702155756374, 0.085488017923844195, 0.02361226606521127, 0.021783497124080092, 0.021166547408477936, 0.020756889979621548, 0.43600177494759196, 0.31028001285571416, 0.18459899380698946, 0.1715284937642364, 0.15199896702758214, 0.12089133548007405, 0.12090406910425358, 0.12104478505889275, 0.12106790756823931, 0.12110204148298398],
    [1.1482483686516933, 0.82519566798847754, 0.60631965941413657, 0.51446949634450378, 0.14587394403416559, 0.05950365694191001, 0.059134109040000285, 0.058952758863629759, 0.058851510175906688, 0.058025644276221568, 0.63439519204264572, 0.33969277991431812, 0.15076883943802938, 0.14833805355329982, 0.14150883379752477, 0.11665758340258567, 0.11670044861110981, 0.11714038215342978, 0.11797150327092823, 0.11782244607808127]]
    """
    list = [[1.1911631306074264, 1.0326002460047061, 0.53642614460666382, 0.34872900130490514, 0.23923569933999683, 0.085452246058788067, 0.02352726964453531, 0.021738678603195331, 0.021092996163994192, 0.020686948031837231, 8.1950040930525527e-07, 8.1903580242208041e-07, 8.1988001991402057e-07, 8.1983260142790491e-07, 8.1975706908470317e-07, 8.1977638825287997e-07, 8.1977378555425437e-07, 8.2014076121620002e-07, 8.1971306468270614e-07, 8.1989344824083702e-07],
    [1.148436774080271, 0.82525790799275955, 0.60642425437999015, 0.51460037018482585, 0.14612762130369858, 0.059859492137331438, 0.059498884211878232, 0.05929892676966949, 0.059194593591231345, 0.05837589416409187, 0.015298565436690512, 0.0020691282391396481, 0.0020802377236021621, 0.0020446736370413805, 0.0020344109706333653, 0.0019040405624723329, 0.0020072422389232313, 0.0020077634325290851, 0.0020079534238409413, 0.0020063415282978924],
    [1.1905957255850586, 1.0319376772851205, 0.53609884278159536, 0.34794440616211608, 0.23835655681233597, 0.084365157214517816, 0.022605139813189373, 0.020782612969124863, 0.02003254891012244, 0.019629620917951472, 0.010467049025574722, 0.0063884967401175458, 0.0066210159578504041, 0.0070210843830951705, 0.0069866274188263881, 0.007045401685077975, 0.0070396672843571883, 0.0070722214878722513, 0.0070641067631039296, 0.0070720851372569919],
    [1.1482911663257593, 0.82513329491636678, 0.60635777019730408, 0.51455836598621829, 0.14588234219916402, 0.059528894179049871, 0.059118799261423369, 0.059016439134619583, 0.058913516245597874, 0.058067972957234613, 0.054618077103265661, 0.028409191108999511, 0.010755596666589578, 0.010970715504907625, 0.011106097370857757, 0.010999406380117874, 0.01095208100603795, 0.010965737657346559, 0.010986112709939056, 0.01098870860213497],
    [1.1913283421824037, 1.0327559077760005, 0.53654019028282973, 0.34868300082386783, 0.23963101953447249, 0.085674919160770219, 0.023761644843105074, 0.021942349693405232, 0.021290676508428516, 0.020914049641599677, 0.033069451717722263, 0.010754743627248181, 0.010787337254451116, 0.011554704925096509, 0.011572036061203893, 0.012031313116868104, 0.012032317605502826, 0.012048897208385612, 0.012054306373108515, 0.012050568802416325],
    [1.1483072925042086, 0.82528262420649523, 0.60646647321297353, 0.51465303945382146, 0.14591990259362678, 0.059660082772301469, 0.059252556672470744, 0.059104513149683714, 0.059002400052374414, 0.058175747419113244, 0.16440477694175165, 0.082807794641487159, 0.025411899716723911, 0.02571224717616729, 0.027093758558074034, 0.026860440373967242, 0.027200482456062852, 0.027219949294459804, 0.027294576105001127, 0.027316822115640987],
    [1.1917253552481681, 1.0332458567323759, 0.53679005646251576, 0.34895927315309427, 0.24031085268105157, 0.086477975007239374, 0.024462517753516922, 0.02265971088192667, 0.022079270563898169, 0.021810437832036511, 0.057457537212017074, 0.024211712514207874, 0.019509624121337499, 0.020871640409407203, 0.020932759573592356, 0.021613516016122887, 0.021629520309736559, 0.021671301435961854, 0.021658207424847913, 0.021772225411051067],
    [1.1481774939317586, 0.82553884546472656, 0.60629425769467571, 0.51472072675389668, 0.14593745107211167, 0.059545729853591564, 0.059172270926414104, 0.059015056869903867, 0.058910279585077621, 0.058109338966537971, 0.28993185952148476, 0.11473524833237217, 0.030290257229892693, 0.030983738898988054, 0.031930471651107364, 0.032277217034677556, 0.033014730477768817, 0.033134037031715365, 0.033210672949525094, 0.03329938383298147],
    [1.191237701041054, 1.0327268907854805, 0.53632485912257266, 0.3484875057512507, 0.23940303508828875, 0.085402852469949941, 0.023641827703683799, 0.021790153178631866, 0.021126183701824945, 0.020805137828906767, 0.084503563444449564, 0.028895138826257862, 0.017959828544788333, 0.020246737920171985, 0.020375669421820071, 0.021118099492220294, 0.020909185894339977, 0.02125558942258501, 0.021300375393050054, 0.021581612735136979],
    [1.1483502839043065, 0.82526150772582818, 0.60638857072423646, 0.51453161172679607, 0.14600387362840994, 0.059730734464464329, 0.059378154514995597, 0.059159805598723657, 0.059057494344170361, 0.058229460783229256, 0.37346286256472488, 0.16161110976560072, 0.02903035760976717, 0.030486840122779517, 0.031484067524785997, 0.033469661525786884, 0.035022919810559744, 0.035400395685706888, 0.035423200965323157, 0.035555433301871525]]
    """
    for k in range(len(list)):
        if k%2 == 0: wpUE.append(list[k][numexp/2:])
        else: wpSO.append(list[k][numexp/2:])
    odUE = list[0][:numexp/2]
    odSO = list[1][:numexp/2]
    index = [10*i for i in range(1,11)]
    
    plt.plot(index, odUE, '-or', label='With OD flows')
    color = ['m', 'c', 'b', 'k', 'g']
    num_wps = [80, 40, 20, 10, 5]
    for i in range(5):
        plt.plot(index, wpUE[i], '-o'+color[i], label='With {} cells'.format(num_wps[i]))
    plt.title('Path flow errors for network in UE')
    plt.xlabel('Percentage of links observed (%)')
    plt.ylabel('Relative error')
    plt.yscale('log')
    plt.legend(loc=0)
    plt.show()
    
    plt.plot(index, odSO, '-or', label='With OD flows')
    for i in range(5):
        plt.plot(index, wpSO[i], '-o'+color[i], label='With {} cells'.format(num_wps[i]))
    plt.title('Path flow errors for network in SO')
    plt.xlabel('Percentage of links observed (%)')
    plt.ylabel('Relative error')
    plt.yscale('log')
    plt.legend(loc=0)
    plt.show()
    
    color = ['c', 'k', 'g', 'b', 'r']
    labels = ['5', '10', '20', '40', '80']
    index2 = range(5)
    for i in range(10):
        if (i+1)%2==0:
            plt.plot(index2, [odUE[i]]*5, '-'+color[i/2], label='{}% of links'.format(index[i]))
            plt.plot(index2, [wpUE[-(k+1)][i] for k in range(5)], '-o'+color[i/2])
    plt.title('Path flow errors for network in UE')
    plt.xlabel('Number of waypoints')
    plt.ylabel('Relative error')
    plt.yscale('log')
    plt.xticks(index2, labels)
    plt.legend(loc=0)
    plt.show()
    
    for i in range(10):
        if (i+1)%2==0:
            plt.plot(index2, [odSO[i]]*5, '-'+color[i/2], label='{}% of links'.format(index[i]))
            plt.plot(index2, [wpSO[-(k+1)][i] for k in range(5)], '-o'+color[i/2])
    plt.title('Path flow errors for network in SO')
    plt.xlabel('Number of waypoints')
    plt.ylabel('Relative error')
    plt.yscale('log')
    plt.xticks(index2, labels)
    plt.legend(loc=0)
    plt.show()


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
    

def display_results_2():
    # FORMAT of output file
    # The first 'N=10' lines are link flow error
    # The second 'N=10' lines are route flow error
    # With each set of N lines, we have 2*len(data) lines, alternating between
    # results without ODs and with ODs
    D = len(data)

    # collect results for UE-type behavior
    results = open('ISTTT_results/ISTTT_results_UE.txt', 'r').readlines()[10:]
    results  = [[float(e) for e in r[1:-2].split(', ')] for r in results]
    est_UE_wp_woODs = [results[2*i][10:] for i in range(D)] #
    est_UE_wp_wiODs = [results[2*i+1][10:] for i in range(D)]
    est_UE_lf = [results[2*i][:10] for i in range(D)]
    est_UE_lf_2 = [results[2*i+1][:10] for i in range(D)]
    
    index = [10*i for i in range(1,11)]
    color = ['m', 'c', 'b', 'k', 'g']
    # num_wps = [120, 80, 40, 20, 10]
    num_wps = [d[0] + d[1] + d[3][0][1] for d in data]

    plt.plot(index, est_UE_lf[0], '-or', label='With OD flows')
    for i in range(D):
        plt.plot(index, est_UE_wp_woODs[i], '-o'+color[i],
                 label='With {} cells'.format(num_wps[i]))
    plt.title('Path flow errors for network in UE: OD vs cellpath')
    plt.xlabel('Percentage of links observed (%)')
    plt.ylabel('Relative error')
    plt.yscale('log')
    plt.legend(loc=0)
    plt.show()
    
    plt.plot(index, est_UE_lf[0], '-or', label='With OD flows')
    for i in range(D):
        plt.plot(index, est_UE_wp_wiODs[i], '-o'+color[i],
                 label='With {} cells'.format(num_wps[i]))
    plt.title('Path flow errors for network in UE: OD vs OD+cellpath')
    plt.xlabel('Percentage of links observed (%)')
    plt.ylabel('Relative error')
    plt.yscale('log')
    plt.legend(loc=0)
    plt.show()
    
    
    # collect results for SO-type behavior
    results = open('ISTTT_results/ISTTT_results_SO.txt', 'r').readlines()[10:]
    results  = [[float(e) for e in r[1:-2].split(', ')] for r in results]
    est_SO_wp_woODs = [results[2*i][10:] for i in range(D)] #
    est_SO_wp_wiODs = [results[2*i+1][10:] for i in range(D)]
    est_SO_lf = [results[2*i][:10] for i in range(D)]
    est_SO_lf_2 = [results[2*i+1][:10] for i in range(D)]
    
    plt.plot(index, est_SO_lf[0], '-or', label='With OD flows')
    for i in range(D):
        plt.plot(index, est_SO_wp_woODs[i], '-o'+color[i], label='With {} cells'.format(num_wps[i]))
    plt.title('Path flow errors for network in SO: OD vs cellpath')
    plt.xlabel('Percentage of links observed (%)')
    plt.ylabel('Relative error')
    plt.yscale('log')
    plt.legend(loc=0)
    plt.show()
    
    plt.plot(index, est_SO_lf[0], '-or', label='With OD flows')
    for i in range(D):
        plt.plot(index, est_SO_wp_wiODs[i], '-o'+color[i], label='With {} cells'.format(num_wps[i]))
    plt.title('Path flow errors for network in SO: OD vs OD+cellpath')
    plt.xlabel('Percentage of links observed (%)')
    plt.ylabel('Relative error')
    plt.yscale('log')
    plt.legend(loc=0)
    plt.show()

def main():
    myseed=29347293
    np.random.seed(myseed)
    random.seed(myseed)

    trials = 100
    #synthetic_data()
    #experiment()
    #ratio_wptrajs_usedpaths()

    # ISTTT experiments:
    run_experiments(SO=False, trials=trials)
    run_experiments(SO=True, trials=trials)
    display_results_2()

    #run_QP_ranks(False)
    #display_ranks()
    #display_ratios()

if __name__ == '__main__':
    main()