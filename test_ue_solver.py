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
    g1, g2, g3, g4 = los_angeles(theta, 'Polynomial')
    linkflows1 = ue.solver(g1)
    linkflows2 = ue.solver(g2)
    linkflows3 = ue.solver(g3)
    linkflows4 = ue.solver(g4)
    d.draw_delays(g1)
    d.draw_delays(g2)
    d.draw_delays(g3)
    d.draw_delays(g4)
    

def main():
    #affine()
    #polynomial()
    los_angeles_ue()
    

if __name__ == '__main__':
    main()