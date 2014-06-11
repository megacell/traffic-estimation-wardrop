'''
Created on Apr 25, 2014

@author: jeromethai
'''

from test_graph import small_grid, small_example, los_angeles_map
import draw_graph as d
import matplotlib.pyplot as plt

od_flows1 = [3.0, 3.0, 1.0, 1.0];

def small_graph():
    #graph = small_example()
    graph = small_grid(od_flows1)
    d.draw(graph)
    d.draw_path(graph, (3,4,1))
    #d.draw_path(graph, (1,5,2))
    plt.figure()

def draw_los_angeles():
    graph = los_angeles_map()
    d.draw(graph)
    plt.figure()

def main():
    #small_graph()
    draw_los_angeles()

if __name__ == '__main__':
    main()
