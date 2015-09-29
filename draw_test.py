'''
Created on Apr 25, 2014

@author: jeromethai
'''

from generate_graph import small_example, los_angeles, los_angeles_2
import draw_graph as d
import matplotlib.pyplot as plt
from cvxopt import matrix


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
    #d.draw(graph, nodes=False)
    #print graph.indlinks
    #for k in range(4): d.draw(graph, link_ids[k], nodes=False)
    #d.draw_ODs(graph, 22)
    # link ids for CDC 2015: [65, 86, 87, 88, 99]
    #link_ids = [(6, 5, 1), (36, 16, 1), (23, 22, 1), (16, 15, 1), (17, 36, 1)]
    link_ids = [(43,24,1), (24,37,1), (37,23,1), (23,22,1), (42,17,1), (17,36,1), (36,16,1), (16,15,1), (9,41,1), (41,8,1), (8,35,1), (35,7,1)]
    for id in link_ids: print graph.indlinks[id]
    d.draw(graph, link_ids, nodes=False)
    

def draw_los_angeles_2():
    theta = matrix([0.0, 0.0, 0.0, 0.15])
    graph = los_angeles_2(theta, 'Polynomial')
    graph.visualize()
    link_ids = [key for key,link in graph.links.items() if link.delayfunc.slope>1.0]
    d.draw(graph, link_ids, nodes=False)


def main():
    #small_graph()
    draw_los_angeles()
    #draw_los_angeles_2()
    

if __name__ == '__main__':
    main()
