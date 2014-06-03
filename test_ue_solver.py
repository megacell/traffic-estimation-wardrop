'''
Created on Apr 21, 2014

@author: jeromethai
'''

import ue_solver as ue
from test_graph import small_grid, small_grid2, small_example

def main():
    graph = small_grid2()
    #graph = small_example()
    linkflows = ue.solver(graph)
    print 'links\' indices: ', graph.indlinks
    print 'UE flow: '
    print linkflows
    unusedpaths = ue.unused_paths(graph)
    print 'Unused paths: ', unusedpaths
    ue.save_mat('../Dropbox/Mega_Cell/data/', 'ue_data2', graph)
    graph.visualize(paths=True)
    return graph, linkflows, unusedpaths

if __name__ == '__main__':
    main()