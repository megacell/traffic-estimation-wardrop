'''
Created on Apr 27, 2014

@author: jeromethai
'''

import networkx as nx
import matplotlib.pyplot as plt


def create_networkx_graph(graph):
    G=nx.DiGraph(indedges={})
    G.add_nodes_from(range(1,graph.numnodes))
    G.add_edges_from([(key[0],key[1]) for key in graph.links.keys()])
    i = 0
    for edge in G.edges(): G.graph['indedges'][edge] = i; i+=1
    return G


def draw(graph, G=None):
    if G is None: G = create_networkx_graph(graph)
    nx.draw(G, pos=graph.nodes_position)
    plt.show()
    
    
def draw_path(graph, pathid, G=None):
    if G is None: G = create_networkx_graph(graph)
    edge_color = ['k'] * G.number_of_edges()
    for link in graph.paths[pathid].links: edge_color[G.graph['indedges'][(link.startnode, link.endnode)]] = 'r'
    nx.draw(G, pos=graph.nodes_position, edge_color=edge_color)
    plt.show()

    