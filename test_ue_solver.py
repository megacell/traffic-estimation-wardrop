'''
Created on Apr 21, 2014

@author: jeromethai
'''

import ue_solver as ue
import graph as g
import test_graph as testg

if __name__ == '__main__':
    grid = testg.smallGrid()
    keys, sol = ue.solver(grid, update=True)
    print 'links\' keys: ', keys
    print 'UE flow: '
    print sol
    
    g.visualize(grid, paths=True)