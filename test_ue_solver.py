'''
Created on Apr 21, 2014

@author: jeromethai
'''

import numpy as np
import ue_solver as ue
import draw_graph as d
from test_graph import small_grid, small_example, los_angeles
from cvxopt import matrix

od_flows1 = [3.0, 3.0, 1.0, 1.0];
od_flows2 = [1.0, 1.0, 1.0, 4.0];

def affine():
    graph = small_grid(od_flows2)
    #graph = small_example()
    linkflows_obs, indlinks_obs = [1.46, 3.54], [(2,1,1), (5,4,1)]
    linkflows = ue.solver(graph, linkflows_obs=linkflows_obs, indlinks_obs=indlinks_obs, soft=100.0)
    print 'links\' indices: ', graph.indlinks
    print 'UE flow: '
    print linkflows
    unusedpaths = ue.unused_paths(graph)
    print
    print 'Unused paths: ', unusedpaths
    #ue.save_mat('../Dropbox/Mega_Cell/data/', 'ue_data2', graph)
    graph.visualize(paths=True)
    return graph, linkflows, unusedpaths


def polynomial():
    graph = small_grid(od_flows2, 'Polynomial', [3.0, 2.0, 1.0])
    linkflows_obs, indlinks_obs = [1.617, 3.383], [(2,1,1), (5,4,1)]
    linkflows = ue.solver(graph)
    #linkflows = ue.solver(graph, linkflows_obs=linkflows_obs, indlinks_obs=indlinks_obs)
    print 'UE flow: '
    print linkflows
    unusedpaths = ue.unused_paths(graph)
    print
    print 'Unused paths: ', unusedpaths
    graph.visualize(links=True, paths=True)
    return graph, linkflows, unusedpaths
    
    
def los_angeles_ue():
    theta = matrix([0.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    theta /= np.sum(theta)
    theta *= 0.15
    graph1, graph2, graph3 = los_angeles(theta, 'Polynomial', True)
    linkflows1 = ue.solver(graph1)
    linkflows2 = ue.solver(graph2)
    linkflows3 = ue.solver(graph3)
    #print 'UE flow: '
    #print 'links\' indices: ', graph.indlinks
    #print linkflows
    #graph.visualize(links=True, only_pos_flows=True)
    d.draw_delays(graph1)
    d.draw_delays(graph2)
    d.draw_delays(graph3)
    

def main():
    #affine()
    #polynomial()
    los_angeles_ue()
    

if __name__ == '__main__':
    main()