'''
Created on Apr 21, 2014

@author: jeromethai
'''

import ue_solver as ue
import graph as g
import test_graph as testg

def main():
    grid = testg.small_grid()
    indlinks, linkflows = ue.solver(grid, update=True)
    print 'links\' indices: ', indlinks
    print 'UE flow: '
    print linkflows
    
    g.visualize(grid, paths=True)

if __name__ == '__main__':
    main()