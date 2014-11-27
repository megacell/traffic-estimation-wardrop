'''
Created on Aug 9, 2014

@author: jeromethai
'''

import Waypoints as w
from generate_graph import los_angeles
from generate_paths import find_UESOpaths
from cvxopt import matrix
import path_solver
import numpy as np
import time
    
    
theta = matrix([0.0, 0.0, 0.0, 0.15])
    
    
def example1():
    """First test to understand the utilization of waypoints"""
    R = w.Rectangle((0.0, 0.0, 10.0, 10.0))
    R.populate(100)
    R.draw_waypoints()
    L = w.Line((0.0, 0.0, 10.0, 10.0))
    L.populate(100)
    L.draw_waypoints()
    point = (5.0, 5.0)
    id = R.closest_to_point(point)
    R.draw_waypoints(wps=[('b',[id],'closest')], ps=[('r',[point],'position')])
    id = L.closest_to_point(point)
    L.draw_waypoints(wps=[('b',[id],'closest')], ps=[('r',[point],'position')])
    polyline = [(5.0, 0.0, 5.0, 5.0), ((5.0, 5.0, 10.0, 5.0))]
    start = time.clock()
    ids = R.closest_to_polyline(polyline, 100)
    print 'Found the closest waypoints to polyline in: ', time.clock() - start
    R.draw_waypoints(wps=[('b',ids,'closest')])
    start = time.clock()
    R.build_partition((6,6), 1.0)
    print 'Build a partition of the space in: ', time.clock() - start
    print R.partition
    start = time.clock()
    ids = R.closest_to_polyline(polyline, 100, True)
    print 'Found the closest waypoints to polyline in: ',time.clock() - start
    R.draw_waypoints(wps=[('b',ids,'closest')])
    
    
def generate_wp(demand=3, data=None, draw=False, voronoi=False, path=None):
    """Generate waypoints following the map of L.A. and draw the map
    
    Parameters:
    ----------
    demand: OD demand
    N0: number of background samples
    N1: number of samples on links
    regions: list of regions, regions[k] = (geometry, N_region)
    res: (n1, n2) s.t. the width is divided into n1 cells and the height into n2 cells
    margin: margin around each cell
    
    Return value:
    ------------
    graph: Graph of L.A.
    WP: waypoint object with waypoints following the map of L.A.
    """
    if data is None:
        N0, N1, scale, regions, res, margin = 20, 40, 0.2, [((3.5, 0.5, 6.5, 3.0), 20)], (12,6), 2.0
    else:
        N0, N1, scale, regions, res, margin = data
    graph = los_angeles(theta, 'Polynomial',path=path)[demand]
    WP = w.sample_waypoints(graph, N0, N1, scale, regions)
    WP.build_partition(res, margin)
    if draw: WP.draw_waypoints(graph, voronoi=voronoi)
    return graph, WP
    
    
def fast_search(SO=False, data=None, demand=3):
    """Show the efficiency of fast search
    1. Get used paths in UE/SO (see generate_paths module)
    2. generate waypoints following the map of L.A.
    3. Draw the closest waypoints to random path in the graph
    4. Compare against fast search cpu time
    """
    g, WP = generate_wp(demand, data)
    paths = find_UESOpaths(SO)
    for p in paths: g.add_path_from_nodes(p)
    g.visualize(general=True)
    k = np.random.random_integers(0, g.numpaths-1)
    path_id = g.paths.keys()[k]
    start = time.clock()
    ids = WP.closest_to_path(g, path_id, 20)
    print time.clock() - start
    WP.draw_waypoints(g, [('r',ids,'closest')], path_id = path_id)
    start = time.clock()
    ids = WP.closest_to_path(g, path_id, 20, True)
    print time.clock() - start
    WP.draw_waypoints(g, [('r',ids,'closest')], path_id = path_id, voronoi=True)
    
    
def compute_wp_flow(SO=False, demand=3, random=False, data=None, path=None):
    """Generate map of L.A., UE path_flow, waypoint trajectories
    1. Generate map of L.A. and waypoints with generate_wp
    2. Get used paths in UE/SO (see generate_paths module)
    3. generate waypoints following the map of L.A.
    4. For each path, find closest waypoints
    
    Parameters:
    ----------
    SO: if True compute the SO
    demand: OD demand
    random: if True, generate random UE/SO paths
    data: waypoint density (N0, N1, scale, regions, res, margin) inputs of generate_wp
    
    Return value:
    -------------
    g: Graph object of L.A.
    p_flow: vector of path flows
    path_wps: dictionary of paths with >tol flow with wp trajectory associated {path_id: wp_ids}
    wp_trajs: list of waypoint trajectories with paths along this trajectory [(wp_traj, path_list, flow)]
    """
    g, WP = generate_wp(demand, data, path=path)
    paths = find_UESOpaths(SO, path=path)
    for p in paths: g.add_path_from_nodes(p)
    g.visualize(general=True)
    p_flow = path_solver.solver(g, update=True, SO=SO, random=random)
    path_wps, wp_trajs = WP.get_wp_trajs(g, 20, True)
    print len(path_wps), len(wp_trajs)
    return g, p_flow, path_wps, wp_trajs
    #print len(path_wps), path_wps
    #print len(wp_trajs), wp_trajs


def main():
    #example1()
    data = (20, 40, 0.2, [((3.5, 0.5, 6.5, 3.0), 20)], (12,6), 2.0)
    #data = (10, 20, 0.2, [((3.5, 0.5, 6.5, 3.0), 10)], (10,5), 2.0)
    #data = (5, 10, 0.2, [((3.5, 0.5, 6.5, 3.0), 5)], (6,3), 2.0)
    #data = (3, 5, 0.2, [((3.5, 0.5, 6.5, 3.0), 2)], (4,2), 2.0)
    #data = (1, 3, 0.2, [((3.5, 0.5, 6.5, 3.0), 1)], (2,2), 2.0)
    #generate_wp(data=data, draw=True, voronoi=False)
    #fast_search(data=data)
    compute_wp_flow(SO=False, random=True, data=data)
    

if __name__ == '__main__':
    main()