'''
Created on Apr 25, 2014

@author: jeromethai
'''

import networkx as nx
import matplotlib.pyplot as plt
from test_graph import small_grid


def main():
    grid = small_grid()
    m,n = 2,3
    G=nx.grid_2d_graph(m,n)
    mapping = {(i,j): n*i+j+1 for i in range(m) for j in range(n)}
    pos = grid.nodes_position
    print pos
    G = nx.relabel_nodes(G,mapping)
    nx.draw(G,pos)
    plt.show()


if __name__ == '__main__':
    main()
