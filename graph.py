'''
Created on Apr 18, 2014

@author: jeromethai
'''

from cvxopt import matrix

class Graph:
    """A graph containing nodes and links"""
    def __init__(self, nodes={}, links={}, ODs={}, paths={}, numnodes=0, numlinks=0, numODs=0, numpaths=0, indlinks={}, indods={}, indpaths={}):
        self.nodes = nodes
        self.links = links
        self.ODs = ODs
        self.paths = paths
        self.numnodes = numnodes
        self.numlinks = numlinks
        self.numODs = numODs
        self.numpaths = numpaths
        self.indlinks = indlinks # indexation for matrix generations
        self.indods = indods # indexation for matrix generations
        self.indpaths = indpaths # indexation for matrix generations
        
    
    def add_node(self):
        self.numnodes += 1
        self.nodes[self.numnodes] = Node({}, {}, {}, {})
        
          
    def add_link(self, startnode, endnode, route=1, flow=0.0, delay=0.0, ffdelay=0.0, delayfunc=None):
        formatStr = 'ERROR: node {} doesn\'t exist, graph countains {} nodes.'
        if startnode == endnode: print 'ERROR: self-loop not allowed.'; return
        if startnode < 1 or startnode > self.numnodes: print formatStr.format(startnode, self.numnodes); return
        if endnode < 1 or endnode > self.numnodes: print formatStr.format(endnode, self.numnodes); return
        
        if (startnode, endnode, route) in self.links:
            print 'ERROR: link ({},{},{}) already exists.'.format(startnode, endnode, route); return
        else:
            link = Link(startnode, endnode, route, flow, delay, ffdelay, delayfunc, {})
            self.indlinks[(startnode, endnode, route)] = self.numlinks
            self.numlinks += 1
            self.links[(startnode, endnode, route)] = link
            self.nodes[startnode].outlinks[(startnode, endnode, route)] = link
            self.nodes[endnode].inlinks[(startnode, endnode, route)] = link
            if delayfunc is not None:
                link.ffdelay = delayfunc.ffdelay
                link.delay = delayfunc.compute_delay(link.flow)
                    
   
    def add_od(self, origin, destination, flow=0.0):
        formatStr = 'ERROR: node {} doesn\'t exist, graph countains {} nodes.'
        if origin == destination: print 'ERROR: self-loop not allowed.'; return
        if origin < 1 or origin > self.numnodes: print formatStr.format(origin, self.numnodes); return
        if destination < 1 or destination > self.numnodes: print formatStr.format(destination, self.numnodes); return
        
        if (origin, destination) in self.ODs:
            print 'ERROR: OD ({},{}) already exists'.format(origin, destination); return
        else:
            self.indods[(origin, destination)] = self.numODs
            self.numODs += 1
            od = OD(origin, destination, flow, {})
            self.ODs[(origin, destination)] = od
            self.nodes[origin].startODs[(origin, destination)] = od
            self.nodes[destination].endODs[(origin, destination)] = od
   
   
    def add_path(self, link_ids):
        
        origin = link_ids[0][0]
        destination = link_ids[len(link_ids)-1][1]
        if not (origin, destination) in self.ODs: print 'ERROR: OD ({},{}) doesn\'t exist.'.format(origin, destination); return
        
        for i in range(len(link_ids)-1):
            if link_ids[i][1] != link_ids[i+1][0]: print 'ERROR: path not valid.'; return
        
        for path in self.ODs[(origin, destination)].paths.values():
            if [(link.startnode,link.endnode,link.route) for link in path.links] == link_ids: print 'ERROR: path already exists.'; return
        
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
        
        
    def visualize(self, general=False, nodes=False, links=False, ODs=False, paths=False):
    
        if general==True:
            print 'Nodes: ', self.nodes
            print 'Number of nodes: ', self.numnodes
            print 'Links: ', self.links
            print 'Number of links: ', self.numlinks
            print 'OD pairs: ', self.ODs
            print 'Number of OD pairs: ', self.numODs
            print 'Paths: ', self.paths
            print 'Number of paths: ', self.numpaths
            print 'Link indexation', self.indlinks
            print 'OD indexation', self.indods
            print 'Path indexation', self.indpaths
            print
  
        if nodes==True:
            for id, node in self.nodes.items():
                print 'Node id: ', id
                print 'In-links: ', node.inlinks
                print 'Out-links: ', node.outlinks
                print 'Start ODs: ', node.startODs
                print 'End ODs: ', node.endODs
                print
     
        if links==True:    
            for id, link in self.links.items():
                print 'Link id: ', id
                print 'Flow: ', link.flow
                print 'Number of paths: ', link.numpaths
                print 'Paths: ', link.paths
                print 'Delay: ', link.delay
                print 'Free flow delay: ', link.ffdelay
                print
        
        if ODs==True:
            for id, od in self.ODs.items():
                print 'OD pair id: ', id
                print 'Flow: ', od.flow
                print 'Number of paths: ', od.numpaths
                print 'Paths: ', od.paths
                print
     
        if paths==True:   
            for id, path in self.paths.items():
                print 'Path id: ', id
                print 'Links: ', [(link.startnode, link.endnode, link.route) for link in path.links]
                print 'Flow: ', path.flow
                print 'Delay: ', path.delay
                print 'Free flow delay: ', path.ffdelay
                print 
        
        
    def get_linkflows(self):
        linkflows = matrix(0.0, (self.numlinks, 1))
        for id,link in self.links.items(): linkflows[self.indlinks[id]] = link.flow
        return linkflows
        
    
class Link:
    """A link in the graph"""
    def __init__(self, startnode, endnode, route=1, flow=0.0, delay=0.0, ffdelay=0.0, delayfunc=None, paths={}, numpaths=0):
        self.startnode = startnode
        self.endnode = endnode
        self.route = route  #if multiple edges
        self.flow = flow  #flow on the link
        self.delay = delay  #delay on the link
        self.ffdelay = ffdelay #free flow delay
        self.delayfunc = delayfunc
        self.paths = paths  #set of paths passing through
        self.numpaths = numpaths
        

class Node:
    """A node in the graph"""
    def __init__(self, inlinks, outlinks, startODs, endODs):
        self.inlinks = inlinks
        self.outlinks = outlinks
        self.startODs = startODs   #set of OD pairs with origin at node
        self.endODs = endODs   #set of OD pairs with destination at node 
        
        
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
    def __init__(self, origin, destination, flow=0.0, paths={}, numpaths=0):
        self.o = origin
        self.d = destination
        self.flow = flow
        self.paths = paths # set of all paths for the OD pair
        self.numpaths = numpaths
        
        
class AffineDelay:
    """Affine Delay function"""
    def __init__(self, ffdelay, slope, type='Affine'):
        self.ffdelay = ffdelay
        self.slope = slope
        self.type = type
        
    def compute_delay(self, flow):
        return self.ffdelay + self.slope*flow
        
        
def create_grid(m, n, inright=None, indown=None, outright=None, outdown=None, 
         inrdelay=None, inddelay=None, outrdelay=None, outddelay=None):
    """construct a grid with m rows and n columns
                                        1 - 2 - 3 - 4
       for a m=2 X n=4 grid, nodes are: |   |   |   |
                                        5 - 6 - 7 - 8
                                              
    Parameters
    ----------
    
    inright: list [k1, k2, k3, k4, k5, ...]: creates ki in-links of node i connected to the right neighbor of i
    indown: list [k1, k2, k3, k4, k5, ...]: creates ki in-links of node i connected to the down neighbor of i
    outright: list [k1, k2, k3, k4, k5, ...]: creates ki out-links of node i connected to the right neighbor of i
    outdown: list [k1, k2, k3, k4, k5, ...]: creates ki out-links of node i connected to the down neighbor of i
    inrdelay: list [[(f11,a11), (f12,a12), ...], [(f21,a21), (f22,a22), ...], ...]: 

    
    """
    grid = Graph()
    [grid.add_node() for i in range(m) for j in range(n)]
    
    if not inright is None:
        if inrdelay is None:
            [grid.add_link(i*n+j+2, i*n+j+1, k+1) for i in range(m) for j in range(n-1) for k in range(inright[i*n+j])]
        else:                
            [grid.add_link(i*n+j+2, i*n+j+1, k+1, delayfunc=AffineDelay(inrdelay[i*n+j][k][0], inrdelay[i*n+j][k][1])) 
             for i in range(m) for j in range(n-1) for k in range(inright[i*n+j])]
            
    if not indown is None:
        if inddelay is None:
            [grid.add_link((i+1)*n+j+1, i*n+j+1, k+1) for i in range(m-1) for j in range(n) for k in range(indown[i*n+j])]
        else:
            [grid.add_link((i+1)*n+j+1, i*n+j+1, k+1, delayfunc=AffineDelay(inddelay[i*n+j][k][0], inddelay[i*n+j][k][1])) 
             for i in range(m-1) for j in range(n) for k in range(indown[i*n+j])]
        
    if not outright is None:
        if outrdelay is None:
            [grid.add_link(i*n+j+1, i*n+j+2, k+1) for i in range(m) for j in range(n-1) for k in range(outright[i*n+j])]
        else:
            [grid.add_link(i*n+j+1, i*n+j+2, k+1, delayfunc=AffineDelay(outrdelay[i*n+j][k][0], outrdelay[i*n+j][k][1]))
             for i in range(m) for j in range(n-1) for k in range(outright[i*n+j])]
        
    if not outdown is None:
        if outddelay is None:
            [grid.add_link(i*n+j+1, (i+1)*n+j+1, k+1) for i in range(m-1) for j in range(n) for k in range(outdown[i*n+j])]
        else:
            [grid.add_link(i*n+j+1, (i+1)*n+j+1, k+1, delayfunc=AffineDelay(outddelay[i*n+j][k][0], outddelay[i*n+j][k][1])) 
             for i in range(m-1) for j in range(n) for k in range(outdown[i*n+j])]
        
    return grid


