'''
Created on Apr 22, 2014

@author: jeromethai
'''

def constraints(graph, linkflows, indlinks):
    """Construct constraints for feasible path flows"""
    return indpaths, C, d, A, d, Aeq, beq


def solver(graph, linkflows, indlinks, update=False, type='lls'):
    """Find a feasible path flow"""
    return indpaths, pathflows