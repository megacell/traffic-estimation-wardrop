'''
Created on Apr 18, 2014

@author: jeromethai
'''

from cvxopt import matrix
import numpy as np
import logging

class Graph:
    """Class Graph containing nodes, links, ODs, paths for traffic assignment"""
    def __init__(self, description=None):
        self.description = description
        self.nodes = {}
        self.links = {}
        self.ODs = {}
        self.paths = {}
        self.numnodes = 0
        self.numlinks = 0
        self.numODs = 0
        self.numpaths = 0
        self.nodes_position = {}
        self.indlinks = {} # indexation for matrix generations
        self.indods = {} # indexation for matrix generations
        self.indpaths = {} # indexation for matrix generations
        
    
    def add_node(self, position=None):
        """Add a node with coordinates as a tuple"""
        self.numnodes += 1
        self.nodes_position[self.numnodes] = position
        self.nodes[self.numnodes] = Node(position)
        
        
    def add_nodes_from_list(self, list):
        """Add nodes from list of positions"""
        for position in list: self.add_node(position)
        
          
    def add_link(self, startnode, endnode, route=1, flow=0.0, delay=0.0, ffdelay=0.0, delayfunc=None):
        """Add a link"""
        formatStr = 'ERROR: node {} doesn\'t exist, graph countains {} nodes.'
        startnode, endnode, route = int(startnode), int(endnode), int(route)
        if startnode == endnode: logging.error('self-loop not allowed.'); return
        if startnode < 1 or startnode > self.numnodes: logging.info(formatStr.format(startnode, self.numnodes)); return
        if endnode < 1 or endnode > self.numnodes: logging.info(formatStr.format(endnode, self.numnodes)); return
        
        if (startnode, endnode, route) in self.links:
            logging.error('link ({},{},{}) already exists.'.format(startnode, endnode, route)); return
        else:
            link = Link(startnode, endnode, route, float(flow), float(delay), float(ffdelay), delayfunc)
            self.indlinks[(startnode, endnode, route)] = self.numlinks
            self.numlinks += 1
            self.links[(startnode, endnode, route)] = link
            self.nodes[startnode].outlinks[(startnode, endnode, route)] = link
            self.nodes[endnode].inlinks[(startnode, endnode, route)] = link
            if delayfunc is not None:
                link.ffdelay = delayfunc.ffdelay
                link.delay = delayfunc.compute_delay(link.flow)
                    
   
    def add_links_from_list(self, list, delaytype):
        """Add links from list
        the list is must contain starnode, endnode, ffdelay, and parameters of delay functions
        """
        for startnode, endnode, route, ffdelay, parameters in list:
            self.add_link(startnode, endnode, route, 0.0, ffdelay, ffdelay, delayfunc=create_delayfunc(delaytype, parameters))
   
   
    def add_od(self, origin, destination, flow=0.0):
        """Add an OD pair"""
        formatStr = 'ERROR: node {} doesn\'t exist, graph countains {} nodes.'
        origin, destination = int(origin), int(destination)
        if origin == destination: logging.error('self-loop not allowed.'); return
        if origin < 1 or origin > self.numnodes: logging.info(formatStr.format(origin, self.numnodes)); return
        if destination < 1 or destination > self.numnodes: logging.info(formatStr.format(destination, self.numnodes)); return
        
        if (origin, destination) in self.ODs:
            logging.error('OD ({},{}) already exists'.format(origin, destination)); return
        else:
            self.indods[(origin, destination)] = self.numODs
            self.numODs += 1
            od = OD(origin, destination, float(flow))
            self.ODs[(origin, destination)] = od
            self.nodes[origin].startODs[(origin, destination)] = od
            self.nodes[destination].endODs[(origin, destination)] = od
   
   
    def add_ods_from_list(self, list):
        """Add OD's from list
        the list is must contain origin, destimation, ffdelay, flow
        """
        for origin, destination, flow in list: self.add_od(origin, destination, flow)
   
   
    def add_path(self, link_ids):
        """Add a path with link_ids a list of link ids"""
        origin = link_ids[0][0]
        destination = link_ids[len(link_ids)-1][1]
        if not (origin, destination) in self.ODs: logging.error('OD ({},{}) doesn\'t exist.'.format(origin, destination)); return
        
        for i in range(len(link_ids)-1):
            if link_ids[i][1] != link_ids[i+1][0]: logging.error('path not valid.'); return
        
        for path in self.ODs[(origin, destination)].paths.values():
            if [(link.startnode,link.endnode,link.route) for link in path.links] == link_ids: logging.error('path already exists.'); return
        
        links = []; delay = 0.0; ffdelay = 0.0
        for id in link_ids:
            link = self.links[(id[0], id[1], id[2])]; links.append(link); delay += link.delay; ffdelay += link.ffdelay
            
        self.ODs[(origin, destination)].numpaths += 1
        route = self.ODs[(origin, destination)].numpaths
        path = Path(origin, destination, route, links, 0.0, delay, ffdelay)
        self.indpaths[(origin, destination, route)] = self.numpaths
        self.numpaths += 1
        self.paths[(origin, destination, route)] = path
        self.ODs[(origin, destination)].paths[(origin, destination, route)] = path
        for link in links:
            self.links[(link.startnode, link.endnode, link.route)].numpaths += 1
            self.links[(link.startnode, link.endnode, link.route)].paths[(origin, destination, route)] = path   
        
    
    def add_path_from_nodes(self, node_ids):
        """Add a path from a list of nodes"""
        link_ids = [(node_ids[k], node_ids[k+1], 1) for k in range(len(node_ids)-1)]
        self.add_path(link_ids)
    
        
    def visualize(self, general=True, nodes=False, links=False, ODs=False, paths=False, only_pos_flows=False, tol=1e-3):
        """Visualize graph"""
        if general:
            logging.info('Description: ', self.description)
            #print 'Nodes: ', self.nodes
            logging.info('Number of nodes: ', self.numnodes)
            #print 'Links: ', self.links
            logging.info('Number of links: ', self.numlinks)
            #print 'OD pairs: ', self.ODs
            logging.info('Number of OD pairs: ', self.numODs)
            #print 'Paths: ', self.paths
            logging.info('Number of paths: ', self.numpaths)
            #print 'Nodes\' position: ', self.nodes_position
            #print 'Link indexation', self.indlinks
            #print 'OD indexation', self.indods
            #print 'Path indexation', self.indpaths

        if nodes:
            for id, node in self.nodes.items():
                logging.info('Node id: ', id)
                logging.info('Position', node.position)
                logging.info('In-links: ', node.inlinks)
                logging.info('Out-links: ', node.outlinks)
                logging.info('Start ODs: ', node.startODs)
                logging.info('End ODs: ', node.endODs)
                logging.info()

        if links:    
            for id, link in self.links.items():
                if link.flow > tol or not only_pos_flows:
                    logging.info('Link id: ', id)
                    logging.info('Flow: ', link.flow)
                    logging.info('Number of paths: ', link.numpaths)
                    logging.info('Paths: ', link.paths)
                    logging.info('Delay: ', link.delay)
                    logging.info('Free flow delay: ', link.ffdelay)
                    if link.delayfunc is not None: logging.info('Type of delay function: ', link.delayfunc.type)
                    logging.info()

        if ODs:
            for id, od in self.ODs.items():
                logging.info('OD pair id: ', id)
                logging.info('Flow: ', od.flow)
                logging.info('Number of paths: ', od.numpaths)
                logging.info('Paths: ', od.paths)
                logging.info()

        if paths:   
            for id, path in self.paths.items():
                logging.info('Path id: ', id)
                logging.info('Links: ', [(link.startnode, link.endnode, link.route) for link in path.links])
                logging.info('Flow: ', path.flow)
                logging.info('Delay: ', path.delay)
                logging.info('Free flow delay: ', path.ffdelay)
                logging.info()

        
        
    def get_linkflows(self):
        """Get link flows in a column cvxopt matrix"""
        linkflows = matrix(0.0, (self.numlinks, 1))
        for id,link in self.links.items(): linkflows[self.indlinks[id]] = link.flow
        return linkflows
    
    
    def get_ffdelays(self):
        """Get ffdelays in a column cvxopt matrix"""
        ffdelays = matrix(0.0, (self.numlinks, 1))
        for id,link in self.links.items(): ffdelays[self.indlinks[id]] = link.delayfunc.ffdelay
        return ffdelays
    
    
    def get_slopes(self):
        """Get slopes in a column cvxopt matrix"""
        slopes = matrix(0.0, (self.numlinks, 1))
        for id,link in self.links.items(): slopes[self.indlinks[id]] = link.delayfunc.slope
        return slopes
    
    
    def get_coefs(self):
        """Get coefficients of the polynomial delay functions in cvxopt.matrix
        
        Return value
        ------------
        coefs[i,j] = coef[j] for link i
        """
        type = self.links.values()[0].delayfunc.type
        if type != 'Polynomial': logging.error('Delay functions must be polynomial'); return
        n, degree = self.numlinks, self.links.values()[0].delayfunc.degree
        coefs = matrix(0.0, (n, degree))
        for id,link in self.links.items():
            for j in range(degree): coefs[self.indlinks[id],j] = link.delayfunc.coef[j]
        return coefs
    
    
    def get_ks(self):
        """Get parameters k1, k2 of the hyperbolic delay function in cvxopt.matrix
        
        Return value
        ------------
        k[i,j] = kj for link i with j=1,2
        """
        type = self.links.values()[0].delayfunc.type
        if type != 'Hyperbolic': logging.error('Delay functions must be hyperbolic'); return
        ks = matrix(0.0, (self.numlinks,2))
        for id,link in self.links.items():
            i = self.indlinks[id]
            ks[i,0], ks[i,1] = link.delayfunc.k1, link.delayfunc.k2
        return ks
        
        
    def get_parameters(self):
        """Get parameters of the graph"""
        type = self.links.values()[0].delayfunc.type
        if type == 'Polynomial': return self.get_coefs(), type
        if type == 'Hyperbolic': return self.get_ks(), type
        
    
    def update_linkflows_linkdelays(self, linkflows):
        """Update link flows and link delays in Graph object"""
        for id,link in self.links.items():
            flow = linkflows[self.indlinks[id]]
            link.flow, link.delay = flow, link.delayfunc.compute_delay(flow)
        
        
    def update_pathdelays(self):
        """Update path delays in Graph object"""
        for path in self.paths.values(): path.delay = sum([link.delay for link in path.links])
        
        
    def update_pathflows(self, pathflows):
        """Update path flows in Graph object"""
        for id,path in self.paths.items(): path.flow = pathflows[self.indpaths[id]]
        
        
    def get_interior_links(self):
        """Get the ids of links in the interior of the graph"""
        intnode_ids, intlk_ids = [], []
        for id,node in self.nodes.items():
            if len(node.inlinks.items()) == 1: intnode_ids.append(id)
        for id,link in self.links.items():
            if link.endnode not in intnode_ids and link.startnode not in intnode_ids:
                intlk_ids.append(id)                
        return intlk_ids

class Link:
    """A link in the graph"""
    def __init__(self, startnode, endnode, route, flow=0.0, delay=0.0, ffdelay=0.0, delayfunc=None):
        self.startnode = startnode
        self.endnode = endnode
        self.route = route  #if multiple edges
        self.flow = flow  #flow on the link
        self.delay = delay  #delay on the link
        self.ffdelay = ffdelay #free flow delay
        self.delayfunc = delayfunc
        self.paths = {}  #set of paths passing through
        self.numpaths = 0

    def repr(self):
        return (self.startnode,self.endnode,self.route)
        

class Node:
    """A node in the graph"""
    def __init__(self, position):
        self.position = position
        self.inlinks = {}
        self.outlinks = {}
        self.startODs = {}   #set of OD pairs with origin at node
        self.endODs = {}   #set of OD pairs with destination at node 
        
        
class Path:
    """A path in the graph"""
    def __init__(self, origin, destination, route, links, flow=0.0, delay=0.0, ffdelay=0.0):
        self.o = origin #origin node
        self.d = destination #destination node
        self.route = route #index of path associated to this od pair
        self.links = links #set of all links on the path
        self.flow = flow
        self.delay = delay
        self.ffdelay = ffdelay
        
        
class OD:
    """OD pairs"""
    def __init__(self, origin, destination, flow=0.0):
        self.o = origin
        self.d = destination
        self.flow = flow
        self.paths = {} # set of all paths for the OD pair
        self.numpaths = 0
       
       
class PolyDelay:
    """Polynomial Delay function
    delay(x) = ffdelay + sum_{k>=1} coef[k]*(x)^k
    example: coef[k]=ffdelay*theta[k]*slope^k w/ theta type = [0.0, 0.0, 0.0, 0.15]"""
    def __init__(self, ffdelay, slope, coef):
        self.ffdelay = ffdelay
        self.slope = slope # this can be seen as inverse capacities
        self.coef = coef
        self.degree = len(coef)
        self.type = 'Polynomial'
        
    def compute_delay(self, flow):
        """Compute delay"""
        return self.ffdelay + np.dot(self.coef, np.power(flow, range(1,self.degree+1)))
    
    
class HyperDelay:
    """Hyperbolic Delay function
    delay(x) = ffdelay - k1/k2 + k1/(k2-x)
    k2 is the max capacity on link (1000/veh/lane)
    example: k1=a*ffdelay/slope, k2=b/slope w/ (a,b) type = (3.5, 3)"""
    def __init__(self, ffdelay, slope, k1, k2):
        self.ffdelay = ffdelay
        self.slope = slope
        self.k1 = k1
        self.k2 = k2
        self.type = 'Hyperbolic'
        
    def compute_delay(self, flow):
        """Compute delay"""
        k1, k2, ffdelay = self.k1, self.k2, self.ffdelay
        return ffdelay - k1/k2 + k1/(k2-flow)
        

def create_delayfunc(type, parameters=None):
    """Create a Delay function of a specific type"""
    if type == 'None': return None
    if type == 'Polynomial': return PolyDelay(parameters[0], parameters[1], parameters[2])
    if type == 'Hyperbolic': return HyperDelay(parameters[0], parameters[1], parameters[2], parameters[3])
    if type == 'Other': return Other(parameters[0], parameters[1], parameters[2])


def create_graph_from_list(list_nodes, list_links, delaytype, list_ods=None, description=None):
    """Create a graph from a list of of nodes and links
    """
    graph = Graph(description)
    graph.add_nodes_from_list(list_nodes)
    graph.add_links_from_list(list_links, delaytype)
    if list_ods is not None: graph.add_ods_from_list(list_ods)
    return graph
