'''
Created on Apr 21, 2014

@author: jeromethai
'''

import ue_solver as ue
from test_graph import small_grid, small_example

od_flows1 = [3.0, 3.0, 1.0, 1.0];
od_flows2 = [1.0, 1.0, 1.0, 4.0];

def affine():
    graph = small_grid(od_flows2)
    #graph = small_example()
    linkflows = ue.solver(graph)
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
    

def main():
    #affine()
    return polynomial()
    

if __name__ == '__main__':
    main()