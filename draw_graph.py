'''
Created on Apr 27, 2014

@author: jeromethai
'''

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


def create_networkx_graph(graph):
    G=nx.DiGraph(indedges={})
    G.add_nodes_from(range(1,graph.numnodes))
    G.add_edges_from([(key[0],key[1]) for key in graph.links.keys()])
    i = 0
    for edge in G.edges(): G.graph['indedges'][edge] = i; i+=1
    return G


def draw(graph, pathid=None, G=None, width=7, alpha=0.5, edge_color='r'):
    """Draw graph
    
    Parameters
    ----------
    graph: graph object
    pathid: if provides a path id, display it in the graph
    G: networkx equivalent of the graph object
    width: width of the highlights for the path
    alpha: transparency of the highlights
    edge_color: color of the highlights
    """
    pos=graph.nodes_position
    if G is None: G = create_networkx_graph(graph)
    nx.draw(G, pos)
    if pathid is not None:
        edgelist = [(link.startnode, link.endnode) for link in graph.paths[pathid].links]
        nx.draw_networkx_edges(G, pos, edgelist=edgelist, width=width, 
                               alpha=alpha, edge_color=edge_color, arrows=False)
    plt.show()
    
    
def draw_delays(graph, G=None, width=7, alpha=0.5, levels=[1.5, 2.0, 3.0], tol=1e-8):
    """Draw graph with delays
    
    Parameters
    ----------
    graph: graph object
    G = networkx equivalent of the graph object
    width: width of the highlights for the path
    alpha: transparency of the highlights
    levels: 3 levels of intensities
    """
    pos=graph.nodes_position
    if G is None: G = create_networkx_graph(graph)
    nx.draw(G, pos)
    edgelists, colors = [[],[],[],[]], ['b', 'g', 'y', 'r']
    u = [f for f in levels]
    u.append(np.inf)
    l = [1.0+tol]
    for f in levels: l.append(f)
    
    for link in graph.links.values():
        delay, ffdelay, startnode, endnode = link.delay, link.ffdelay, link.startnode, link.endnode
        for i in range(4):
            if ffdelay * l[i] <= delay < ffdelay * u[i]: edgelists[i].append((startnode,endnode))
    
    for i in range(4):          
        nx.draw_networkx_edges(G, pos, edgelist=edgelists[i],
                               width=width, alpha=alpha, edge_color=colors[i],arrows=False)
   
    plt.show()