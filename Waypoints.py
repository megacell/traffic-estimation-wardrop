'''
Created on Aug 9, 2014

@author: jeromethai
'''

import numpy as np
import numpy.random as ra
from util import sample_line, sample_box, create_networkx_graph
import networkx as nx    
import matplotlib.pyplot as plt


class Waypoints:
    """Waypoints containing geometry, N waypoints, and a shape"""
    def __init__(self, geo, shape='Shape'):
        self.geometry = geo
        self.shape = shape
        self.N = 0
        self.wp = {}
        
        
class Rectangle(Waypoints):
    """Rectangle containing geo=(x1,y1,x2,y2), N waypoints, and a shape"""
    def __init__(self, geo):
        Waypoints.__init__(self, geo, 'Rectangle')
        
    def populate(self, N, first=1):
        """Uniformly sample N points in rectangle
        with first the first key used in wp"""
        self.N = N
        ps = sample_box(N, self.geometry)
        self.wp = {id: p for id,p in enumerate(ps,first)}
        if self.shape == 'Bounding box': self.N0 = self.N


class BoundingBox(Rectangle):
    """BoundingBox containing geo=(x1,y1,x2,y2), N waypoints, shape, lines, regions"""
    def __init__(self, geo):
        Rectangle.__init__(self, geo)
        self.shape = 'Bounding box'
        self.lines = {}
        self.num_lines = 0
        self.regions = {}
        self.num_regions = 0
        self.N0 = 0 # number of uniform samples in the whole region
        
    def add_rectangle(self, geo, N):
        """Add a rectangular region with N points"""
        r = Rectangle(geo)
        r.populate(N, self.N+1)
        self.num_regions += 1
        self.regions[self.num_regions] = r
        self.N += N
        self.wp = dict(self.wp.items() + r.wp.items())
        
    def add_line(self, geo, N, scale):
        """Add a line with N points"""
        l = Line(geo)
        l.populate(N, self.N+1, scale)
        self.num_lines += 1
        self.lines[self.num_lines] = l
        self.N += N
        self.wp = dict(self.wp.items() + l.wp.items())
                    
        
class Line(Waypoints):
    """Class Line containing geo=(x1,y1,x2,y2) waypoints"""
    def __init__(self, geo):
        Waypoints.__init__(self, geo, 'Line')
        
    def populate(self, N, first, scale):
        """Sample N points along line
        with first the first key used in wp"""
        self.N = N
        ps = sample_line(N, self.geometry, scale)
        self.wp = {id: p for id,p in enumerate(ps,first)}


def sample_waypoints(graph, N0, N1, regions, margin):
    """Sample waypoints on graph
    
    Parameters:
    -----------
    graph: Graph object
    N0: number of background samples
    N1: number of samples on links
    regions: list of regions, regions[k] = (geometry, N_region)
    margin: % size of margin around the graph
    """
    xs = [p[0] for p in graph.nodes_position.values()]
    ys = [p[1] for p in graph.nodes_position.values()]
    min_x, max_x, min_y, max_y = min(xs), max(xs), min(ys), max(ys)
    w, h = max_x-min_x, max_y-min_y
    x1, x2, y1, y2 = min_x - w*margin, max_x + w*margin, min_y - h*margin, max_y + h*margin
    WP = BoundingBox((x1, y1, x2, y2))
    WP.populate(N0)
    total_length, lines = 0, []
    for link in graph.links.values():
        xs, ys = graph.nodes_position[link.startnode]
        xt, yt = graph.nodes_position[link.endnode]
        length = np.linalg.norm([xs-xt, ys-yt])
        total_length += length
        lines.append([(xs,ys,xt,yt), length])
    weights = [line[1]/total_length for line in lines]
    Ns = ra.multinomial(N1, weights, size=1)[0]
    for k,line in enumerate(lines): WP.add_line(line[0], Ns[k], 0.1)
    for r in regions: WP.add_rectangle(r[0], r[1])
    return WP


def draw_waypoints(graph, WP, G=None):
    """Draw graph and waypoints
    """
    pos=graph.nodes_position
    if G is None: G = create_networkx_graph(graph)
    nx.draw_networkx_edges(G, pos, arrows=False)
    if WP.N0 > 0:
        xs = [WP.wp[i+1][0] for i in range(WP.N0)]
        ys = [WP.wp[i+1][1] for i in range(WP.N0)]
        plt.plot(xs, ys, 'co', label='uniform')
    if len(WP.lines) > 0:
        xs = [p[0] for line in WP.lines.values() for p in line.wp.values()]
        ys = [p[1] for line in WP.lines.values() for p in line.wp.values()]
        plt.plot(xs, ys, 'mo', label='lines')
    if len(WP.regions) > 0:
        xs = [p[0] for r in WP.regions.values() for p in r.wp.values()]
        ys = [p[1] for r in WP.regions.values() for p in r.wp.values()]
        plt.plot(xs, ys, 'go', label='regions')
    plt.legend()
    plt.show()
    

if __name__ == '__main__':
    pass