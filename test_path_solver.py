'''
Created on Apr 22, 2014

@author: jeromethai
'''

import test_ue_solver as testue
import path_solver as path
from cvxopt import matrix

def main():
    grid, linkflows = testue.main()
    pathflows = path.solver(grid, linkflows)
    print 'paths\' indices: ', grid.indpaths
    print pathflows



if __name__ == '__main__':
    main()