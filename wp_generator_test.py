'''
Created on Aug 9, 2014

@author: jeromethai
'''

import Waypoints as w
from generate_graph import los_angeles
from generate_paths import find_UESOpaths
from cvxopt import matrix
import path_solver as path
import numpy as np
    
theta = matrix([0.0, 0.0, 0.0, 0.15])
    
    
def example():
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
    ids = R.closest_to_polyline(polyline, 100)
    R.draw_waypoints(wps=[('b',ids,'closest')])
    
    
def generate_wp(demand=3, draw=False):
    """Generate waypoints following the map of L.A. and draw the map"""
    graph = los_angeles(theta, 'Polynomial')[demand]
    regions = [((3.5, 0.5, 6.5, 3.0), 50)]
    WP = w.sample_waypoints(graph, 50, 100, regions, 0.05)
    if draw: WP.draw_waypoints(graph)
    return graph, WP
    
    
def test3(SO=False, demand=3):
    """
    1. Get used paths in UE/SO (see generate_paths module)
    2. generate waypoints following the map of L.A.
    3. Draw the closest waypoints to random path in the graph
    """
    g, WP = generate_wp(demand)
    paths = find_UESOpaths(SO)
    for p in paths: g.add_path_from_nodes(p)
    g.visualize(general=True)
    k = np.random.random_integers(0, g.numpaths-1)
    path_id = g.paths.keys()[k]
    ids = WP.closest_to_path(g, path_id, 20)
    WP.draw_waypoints(g, [('r',ids,'closest')], path_id = path_id)
    p_flow = path.solver(g, update=True, SO=SO)
    wp_flow = WP.generate_wp_flows(g, 20)
    print len(wp_flow), wp_flow
    
    
#def test4(SO, ):
    

def main():
    #example()
    #generate_wp(draw=True)
    test3()

if __name__ == '__main__':
    main()