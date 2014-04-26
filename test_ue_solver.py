'''
Created on Apr 21, 2014

@author: jeromethai
'''

import ue_solver as ue
from test_graph import small_grid

def main():
    grid = small_grid()
    linkflows = ue.solver(grid)
    print 'links\' indices: ', grid.indlinks
    print 'UE flow: '
    print linkflows
    unusedpaths = ue.unused_paths(grid)
    print 'Unused paths: ', unusedpaths
    return grid, linkflows, unusedpaths

if __name__ == '__main__':
    main()