'''
Created on Apr 25, 2014

@author: jeromethai
'''

from test_graph import small_example, los_angeles
import draw_graph as d
import matplotlib.pyplot as plt

od_flows1 = [3.0, 3.0, 1.0, 1.0];

def small_graph():
    #graph = small_example()
    graph = small_grid(od_flows1)
    d.draw(graph)
    d.draw(graph, (3,4,1))


def draw_los_angeles():
    graph = los_angeles()[0]
    link_ids = [(17,24,1),(24,40,1),(14,21,1),(16,23,1)]
    #link_ids = [(10,9,1),(19,18,1),(4,5,1),(29,21,1)]
    #link_ids = graph.links.keys()
    d.draw(graph, link_ids)


def main():
    #small_graph()
    draw_los_angeles()

if __name__ == '__main__':
    main()
