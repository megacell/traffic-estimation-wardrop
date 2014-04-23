'''
Created on Apr 22, 2014

@author: jeromethai
'''

import test_ue_solver as testue
import path_solver as path
from cvxopt import matrix

def main():
    grid, indlinks, linkflows = testue.main()
    pathflows = path.solver(grid, linkflows, indlinks)
    print 'links\' indices: ', indlinks
    print 'paths\' indices: ', grid.indpaths
    print pathflows
    #print U
    #print r
    #print C
    #print d


if __name__ == '__main__':
    main()