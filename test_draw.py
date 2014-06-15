'''
Created on Apr 25, 2014

@author: jeromethai
'''

from test_graph import small_grid, small_example, los_angeles
import draw_graph as d
import matplotlib.pyplot as plt

od_flows1 = [3.0, 3.0, 1.0, 1.0];

def small_graph():
    #graph = small_example()
    graph = small_grid(od_flows1)
    d.draw(graph)
    d.draw(graph, (3,4,1))


def draw_los_angeles():
    graph = los_angeles()
    d.draw(graph)


def main():
    small_graph()
    #draw_los_angeles()

if __name__ == '__main__':
    main()
