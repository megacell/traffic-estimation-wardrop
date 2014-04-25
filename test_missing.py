'''
Created on Apr 23, 2014

@author: jeromethai
'''

import path_solver as path
import test_ue_solver as testue
import missing as mis

def main():
    grid, linkflows, unusedpaths = testue.main()
    A = path.incidence(grid)
    ind, misA, misflows = mis.remove_meas(grid, [(2,5,1), (5,2,1)], A, linkflows)
    print ind
    print grid.indlinks
    print misA
    print misflows
    print mis.remove_meas_rand(grid, 2, A, linkflows)

if __name__ == '__main__':
    main()