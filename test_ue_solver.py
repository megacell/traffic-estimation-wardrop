'''
Created on Apr 21, 2014

@author: jeromethai
'''

import ue_solver as ue
import Graph as g
import test_graph as testg

def main():
    grid = testg.small_grid()
    linkflows = ue.solver(grid, update=True)
    print 'links\' indices: ', grid.indlinks
    print 'UE flow: '
    print linkflows
    #g.visualize(grid, links=True, paths=True)
    return grid, linkflows

if __name__ == '__main__':
    main()