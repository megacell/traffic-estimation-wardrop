'''
Created on Apr 22, 2014

@author: jeromethai
'''

import test_ue_solver as testue
import path_solver as path

def main():
    grid, linkflows, unusedpaths = testue.main()
    pathflows = path.solver(grid, linkflows, unusedpaths=unusedpaths)
    print 'links\' indices: ', grid.indlinks
    print 'paths\' indices: ', grid.indpaths
    print pathflows
    print path.vec_feas_paths(grid, unusedpaths=unusedpaths)

if __name__ == '__main__':
    main()