'''
Created on Apr 25, 2014

@author: jeromethai
'''

from test_graph import small_grid, small_example
import draw_graph as d
import matplotlib.pyplot as plt

od_flows1 = [3.0, 3.0, 1.0, 1.0];

def main():
    #graph = small_example()
    graph = small_grid(od_flows1)
    d.draw(graph)
    d.draw_path(graph, (3,4,1))
    #d.draw_path(graph, (1,5,2))
    plt.figure()

if __name__ == '__main__':
    main()
