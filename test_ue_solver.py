'''
Created on Apr 21, 2014

@author: jeromethai
'''

import ue_solver as ue
import graph as g
import test_graph as testg

def main():
    grid = testg.small_grid()
    linkflows = ue.solver(grid, update=True)
    print 'links\' indices: ', grid.indlinks
    print 'UE flow: '
    print linkflows
<<<<<<< HEAD
    #g.visualize(grid, paths=True)
    #g.visualize(grid, links=True, paths=True)
    return grid, linkflows
=======
    
    g.visualize(grid, paths=True)
>>>>>>> parent of 47bc38d... complete path_solver.py and test_path_solver

if __name__ == '__main__':
    main()