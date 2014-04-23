'''
Created on Apr 22, 2014

@author: jeromethai
'''

import test_ue_solver as testue
import path_solver as path
from cvxopt import matrix

def main():
    grid, linkflows = testue.main()
    pathflows = path.solver(grid, linkflows, update=True)
    print 'links\' indices: ', grid.indlinks
    print 'paths\' indices: ', grid.indpaths
    print pathflows
    print path.vec_feas_paths(grid)

if __name__ == '__main__':
    main()