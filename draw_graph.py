'''
Created on Apr 27, 2014

@author: jeromethai
'''

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from util import create_networkx_graph


def draw(graph, link_ids=None, G=None, width=7, alpha=0.5, edge_color='r', nodes=True):
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
    m=graph.numnodes
    pos=graph.nodes_position
    if G is None: G = create_networkx_graph(graph)
    if nodes: 
        nx.draw(G, pos, arrows=False, node_color='w')
    else: 
        nx.draw_networkx_edges(G, pos, arrows=False)
    if link_ids is not None:
        nx.draw_networkx_edges(G, pos, edgelist=[(id[0],id[1]) for id in link_ids], width=width, 
                               alpha=alpha, edge_color=edge_color, arrows=False)
    plt.show()
    

def draw_ODs(graph, dest, G=None):
    """Draw graph and highlights all the ODs that share the same destination
    """
    m=graph.numnodes
    pos=graph.nodes_position
    if G is None: G = create_networkx_graph(graph)
    nx.draw_networkx_edges(G, pos, arrows=False)
    nx.draw_networkx_nodes(G,pos,nodelist=range(1,m+1),node_color='w',node_size=400)
    nx.draw_networkx_nodes(G,pos,nodelist=[dest],node_color='r',node_size=400)
    os = [o for (o,d) in graph.ODs.keys() if d==dest]
    nx.draw_networkx_nodes(G,pos,nodelist=os,node_color='c',node_size=400)
    nx.draw_networkx_labels(G,pos,{k:str(k) for k in range(1,m+1)},font_size=14)
    plt.show()

    
def draw_delays(graph, linkflows=None, G=None, width=7, alpha=0.5, levels=[1.5, 2.0, 3.0], tol=1e-4):
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
    #nx.draw(G, pos, arrows=False, node_color='w')
    nx.draw_networkx_edges(G, pos, arrows=False)
    edgelists, colors = [[],[],[],[]], ['b', 'g', 'y', 'r']
    u = [f for f in levels]
    u.append(np.inf)
    l = [1.0+tol]
    if linkflows is not None:
        l = [1.0]
    else:
        linkflows = [np.inf]*graph.numlinks
    for f in levels: l.append(f)
    #max = 0.0
    for id,link in graph.links.items():
        if linkflows[graph.indlinks[id]] > tol:
            delay, ffdelay, startnode, endnode = link.delay, link.ffdelay, link.startnode, link.endnode
            #if delay/ffdelay > max: max, s, e = delay/ffdelay, startnode, endnode 
            for i in range(4):
                if ffdelay * l[i] <= delay < ffdelay * u[i]: edgelists[i].append((startnode,endnode))
    #print max, s, e
    for i in range(4):          
        nx.draw_networkx_edges(G, pos, edgelist=edgelists[i],
                               width=width, alpha=alpha, edge_color=colors[i],arrows=False)
   
    plt.show()