'''
Created on Apr 18, 2014

@author: jeromethai
'''

from sets import Set
import itertools as iter

class Graph:
    """A graph containing nodes and links"""
    def __init__(self, nodes, links, ODs, paths, numnodes=0, numlinks=0, numODs=0, numpaths=0):
        self.nodes = nodes
        self.links = links
        self.ODs = ODs
        self.paths = paths
        self.numnodes = numnodes
        self.numlinks = numlinks
        self.numODs = numODs
        self.numpaths = numpaths
        
        
    def addNode(self):
        self.numnodes += 1
        self.nodes[self.numnodes] = Node({}, {}, {}, {})
        
          
    def addLink(self, startnode, endnode, route=1, flow=0.0, delay=0.0, ffdelay=0.0, delayfunc=None):
        formatStr = 'ERROR: node {} doesn\'t exist, graph countains {} nodes'
        if startnode == endnode: print 'ERROR: self-loop not allowed'; return
        if startnode < 1 or startnode > self.numnodes: print formatStr.format(startnode, self.numnodes); return
        if endnode < 1 or endnode > self.numnodes: print formatStr.format(endnode, self.numnodes); return
        
        if (startnode, endnode, route) in self.links:
            print 'ERROR: link ({},{},{}) already exists'.format(startnode, endnode, route); return
        else:
            link = Link(startnode, endnode, route, flow, delay, ffdelay, delayfunc, {})
            self.numlinks += 1
            self.links[(startnode, endnode, route)] = link
            self.nodes[startnode].outlinks[(startnode, endnode, route)] = link
            self.nodes[endnode].inlinks[(startnode, endnode, route)] = link
            if not delayfunc is None:
                link.ffdelay = delayfunc.ffdelay
                link.delay = delayfunc.computeDelay(link.flow)
                    
   
    def addOD(self, origin, destination, flow=0.0):
        formatStr = 'ERROR: node {} doesn\'t exist, graph countains {} nodes'
        if origin == destination: print 'ERROR: self-loop not allowed'; return
        if origin < 1 or origin > self.numnodes: print formatStr.format(origin, self.numnodes); return
        if destination < 1 or destination > self.numnodes: print formatStr.format(destination, self.numnodes); return
        
        if (origin, destination) in self.ODs:
            print 'ERROR: OD ({},{}) already exists'.format(origin, destination); return
        else:
            self.numODs += 1
            od = OD(origin, destination, flow, {})
            self.ODs[(origin, destination)] = od
            self.nodes[origin].startODs[(origin, destination)] = od
            self.nodes[destination].endODs[(origin, destination)] = od
   
   
    def addPath(self, linkIds):
        
        origin = linkIds[0][0]
        destination = linkIds[len(linkIds)-1][1]
        if not (origin, destination) in self.ODs: print 'ERROR: OD ({},{}) doesn\'t exist'.format(origin, destination); return
        
        for i in range(len(linkIds)-1):
            if linkIds[i][1] != linkIds[i+1][0]: print 'ERROR: path not valid'; return
        
        for path in self.ODs[(origin, destination)].paths.values():
            if [(link.startnode,link.endnode,link.route) for link in path.links] == linkIds: print 'ERROR: path already exists'; return
        
        links = []; delay = 0.0; ffdelay = 0.0
        for id in linkIds:
            link = self.links[(id[0], id[1], id[2])]; links.append(link); delay += link.delay; ffdelay += link.ffdelay
            
        path = Path(origin, destination, links, 0.0, delay, ffdelay)
        self.numpaths += 1
        self.ODs[(origin, destination)].numpaths += 1
        route = self.ODs[(origin, destination)].numpaths
        self.paths[(origin, destination, route)] = path
        self.ODs[(origin, destination)].paths[(origin, destination, route)] = path
        for link in links:
            self.links[(link.startnode, link.endnode, link.route)].numpaths += 1
            self.links[(link.startnode, link.endnode, link.route)].paths[(origin, destination, route)] = path   
        
    
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
    def __init__(self, origin, destination, links, flow=0.0, delay=0.0, ffdelay=0.0):
        self.o = origin #origin node
        self.d = destination #destination node
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
        
    def computeDelay(self, flow):
        return self.ffdelay + self.slope*flow
        
        
def grid(n, m, inRight=None, inDown=None, outRight=None, outDown=None, 
         inRdelay=None, inDdelay=None, outRdelay=None, outDdelay=None):
    """construct a grid with n rows and m columns"""
    grid = Graph({},{},{},{})
    [grid.addNode() for i in range(n) for j in range(m)]
    
    if not inRight is None:
        if inRdelay is None:
            [grid.addLink(i*m+j+2, i*m+j+1, k+1) for i in range(n) for j in range(m-1) for k in range(inRight[i*m+j])]
        else:                
            [grid.addLink(i*m+j+2, i*m+j+1, k+1, delayfunc=AffineDelay(inRdelay[i*m+j][k][0], inRdelay[i*m+j][k][1])) 
             for i in range(n) for j in range(m-1) for k in range(inRight[i*m+j])]
            
    if not inDown is None:
        if inDdelay is None:
            [grid.addLink((i+1)*m+j+1, i*m+j+1, k+1) for i in range(n-1) for j in range(m) for k in range(inDown[i*m+j])]
        else:
            [grid.addLink((i+1)*m+j+1, i*m+j+1, k+1, delayfunc=AffineDelay(inDdelay[i*m+j][k][0], inDdelay[i*m+j][k][1])) 
             for i in range(n-1) for j in range(m) for k in range(inDown[i*m+j])]
        
    if not outRight is None:
        if outRdelay is None:
            [grid.addLink(i*m+j+1, i*m+j+2, k+1) for i in range(n) for j in range(m-1) for k in range(outRight[i*m+j])]
        else:
            [grid.addLink(i*m+j+1, i*m+j+2, k+1, delayfunc=AffineDelay(outRdelay[i*m+j][k][0], outRdelay[i*m+j][k][1]))
             for i in range(n) for j in range(m-1) for k in range(outRight[i*m+j])]
        
    if not outDown is None:
        if outDdelay is None:
            [grid.addLink(i*m+j+1, (i+1)*m+j+1, k+1) for i in range(n-1) for j in range(m) for k in range(outDown[i*m+j])]
        else:
            [grid.addLink(i*m+j+1, (i+1)*m+j+1, k+1, delayfunc=AffineDelay(outDdelay[i*m+j][k][0], outDdelay[i*m+j][k][1])) 
             for i in range(n-1) for j in range(m) for k in range(outDown[i*m+j])]
        
    return grid


def visualize(graph, general=False, nodes=False, links=False, ODs=False, paths=False):
    
    if general==True:
        print 'Nodes: ', graph.nodes
        print 'Number of nodes: ', graph.numnodes
        print 'Links: ', graph.links
        print 'Number of links: ', graph.numlinks
        print 'OD pairs: ', graph.ODs
        print 'Number of OD pairs: ', graph.numODs
        print 'Paths: ', graph.paths
        print 'Number of paths: ', graph.numpaths
        print
  
    if nodes==True:
        for id, node in graph.nodes.items():
            print 'Node id: ', id
            print 'In-links: ', node.inlinks
            print 'Out-links: ', node.outlinks
            print 'Start ODs: ', node.startODs
            print 'End ODs: ', node.endODs
            print
     
    if links==True:    
        for id, link in graph.links.items():
            print 'Link id: ', id
            print 'Number of paths: ', link.numpaths
            print 'Paths: ', link.paths
            print 'Delay: ', link.delay
            print 'Free flow delay: ', link.ffdelay
            print
        
    if ODs==True:
        for id, od in graph.ODs.items():
            print 'OD pair id: ', id
            print 'Flow: ', od.flow
            print 'Number of paths: ', od.numpaths
            print 'Paths: ', od.paths
            print
     
    if paths==True:   
        for id, path in graph.paths.items():
            print 'Path id: ', id
            print 'Origin: ', path.o
            print 'Destination: ', path.d
            print 'Links: ', [(link.startnode, link.endnode, link.route) for link in path.links]
            print 'Delay: ', path.delay
            print 'Free flow delay: ', path.ffdelay
            print 