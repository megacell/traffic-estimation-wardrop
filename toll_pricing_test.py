'''
Created on Sep 9, 2014

@author: jeromethai
'''

import toll_pricing as toll
from generate_graph import los_angeles
from cvxopt import matrix

theta = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
degree = len(theta)
graph = los_angeles(theta, 'Polynomial')[3]


def test():
    ws_so = [1e-8, 1e-6, 1e-4, 1e-2]
    ws_toll = [1e-8, 1e-6, 1e-4, 1e-2]
    costs, toll_collection, UE_cost, SO_cost = toll.main_solver(graph, theta, ws_so, ws_toll)
    print costs-4700.0
    print toll_collection
    print 'UE_cost:', UE_cost
    print 'SO_cost:', SO_cost


def main():
    test()


if __name__ == '__main__':
    main()