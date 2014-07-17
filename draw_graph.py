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


def draw(graph, link_ids=None, G=None, width=7, alpha=0.5, edge_color='r'):
    """Draw graph
    
    Parameters
    ----------
    graph: graph object
    link_ids: if provided, highlight links
    G: networkx equivalent of the graph object
    width: width of the highlights for the path
    alpha: transparency of the highlights
    edge_color: color of the highlights
    """
    pos=graph.nodes_position
    if G is None: G = create_networkx_graph(graph)
    nx.draw(G, pos)
    if link_ids is not None:
        nx.draw_networkx_edges(G, pos, edgelist=[(id[0],id[1]) for id in link_ids], width=width, 
                               alpha=alpha, edge_color=edge_color, arrows=False)
    plt.show()
    
    
def draw_delays(graph, linkflows=None, G=None, width=7, alpha=0.5, levels=[1.5, 2.0, 3.0], tol=1e-8):
    """Draw graph with delays
    
    Parameters
    ----------
    graph: graph object
    linkflows: only color edges with >0 flow
    G: networkx equivalent of the graph object
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
    if linkflows is not None:
        l = [1.0]
    else:
        linkflows = [np.inf]*graph.numlinks
    for f in levels: l.append(f)
    
    for id,link in graph.links.items():
        if linkflows[graph.indlinks[id]] > tol:
            delay, ffdelay, startnode, endnode = link.delay, link.ffdelay, link.startnode, link.endnode
            for i in range(4):
                if ffdelay * l[i] <= delay < ffdelay * u[i]: edgelists[i].append((startnode,endnode))
    
    for i in range(4):          
        nx.draw_networkx_edges(G, pos, edgelist=edgelists[i],
                               width=width, alpha=alpha, edge_color=colors[i],arrows=False)
   
    plt.show()