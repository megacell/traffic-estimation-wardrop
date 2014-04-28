'''
Created on Apr 22, 2014

@author: jeromethai
'''

import test_ue_solver as testue
import path_solver as path

def main():
    graph, linkflows, unusedpaths = testue.main()
    pathflows = path.solver(graph, linkflows, unusedpaths=unusedpaths)
    print 'links\' indices: ', graph.indlinks
    print 'paths\' indices: ', graph.indpaths
    print pathflows
    print path.vec_feas_paths(graph, unusedpaths=unusedpaths)

if __name__ == '__main__':
    main()