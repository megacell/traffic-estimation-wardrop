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
    link_ids = {}
    link_ids[0] = graph.indlinks.keys()
    link_ids[1] = [(36,37,1), (13,14,1), (17,8,1), (24,17,1), (28,22,1), (14,13,1), (17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    link_ids[2] = [(17,24,1), (24,40,1), (14,21,1), (16,23,1)]
    link_ids[3] = [(10,9,1), (19,18,1), (4,5,1), (29,21,1)]
    #d.draw(graph)
    #d.draw(graph, link_ids)
    #d.draw(graph, nodes=False)
    for k in range(4): d.draw(graph, link_ids[k], nodes=False)
    #d.draw_ODs(graph, 22)


def main():
    #small_graph()
    draw_los_angeles()

if __name__ == '__main__':
    main()
