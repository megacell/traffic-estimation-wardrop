'''
Created on Apr 22, 2014

@author: jeromethai
'''
def remove_meas(graph, link_ids, indlinks, A, linkflows):
    """Remove measurements from links with ids link_ids
    
    Parameters
    ----------
    graph: graph object
    link_ids: list of ids [id1 , id2, id3, ...] with ids of the form id = (link.startnode, link.endnode, link.route)
    indlinks: dictionary of link id to index
    A: matrix of incidence links-paths
    linkflows: matrix of link flows
    """
    ind = range(graph.numlinks)
    for linkid in link_ids: ind.remove(indlinks[linkid])
    return ind, A[ind,:], linkflows[ind,:]


def remove_random_meas(graph, k, indlinks, A, linkflows):
    """Remove randomly k measurements from links"""
    return