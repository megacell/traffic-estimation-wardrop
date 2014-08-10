'''
Created on Aug 9, 2014

@author: jeromethai
'''

import Waypoints as w
from generate_graph import los_angeles

    
def test1():
    """First test to understand the utilization of waypoints"""
    bbox = w.BoundingBox((0.0, 0.0, 1.0, 1.0))
    bbox.populate(3)
    print bbox.N
    print bbox.wp
    bbox.add_rectangle((0.0, 0.0, 0.1, 0.1), 3)
    print bbox.N
    print bbox.wp
    print bbox.regions
    print bbox.regions.values()[0].wp
    print bbox.regions.values()[0].N
    bbox.add_line((0.5, 0.5, 1.0, 1.0), 3, 0.1)
    print bbox.N
    print bbox.wp
    print bbox.lines
    print bbox.lines.values()[0].wp
    print bbox.lines.values()[0].N
    print bbox.N0
    
    
def test2():
    """Generate waypoints following the map of L.A. and draw the map"""
    graph = los_angeles()[0]
    WP = w.sample_waypoints(graph, 50, 100, [((3.5, 0.5, 6.5, 3.0), 50)], 0.05)
    print len(WP.wp)
    w.draw_waypoints(graph, WP)
    
    

def main():
    #test1()
    test2()

if __name__ == '__main__':
    main()