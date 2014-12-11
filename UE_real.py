# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 16:15:19 2014

@author: hugo
"""

import sys
sys.path.append("/home/hugo/Desktop/Hugo/UE_Jerome")
import numpy as np
import ue_solver as ue
import draw_graph as d
from generate_graph import los_angeles_2
from cvxopt import matrix, mul
import time

def test2(delaytype):
    if delaytype == 'Polynomial': theta = matrix([0.0, 0.0, 0.0, 0.15, 0.0, 0.0])
    if delaytype == 'Hyperbolic': theta = (3.5, 3.0)
    g = los_angeles_2(theta, delaytype)
    l, x = ue.solver(g, update=True, full=True)
    d.draw_delays(g)
    print (max(mul(l,g.get_slopes())))
    print ('cost UE:', sum([link.delay*link.flow for link in g.links.values()]))
    l2, x2 = ue.solver(g, update=True, full=True, SO=True)
    print (max(mul(l2,g.get_slopes())))
    print ('cost SO:', sum([link.delay*link.flow for link in g.links.values()]))

def test_osm_box(delaytype):
    t0=time.time()
    box = [34.124918, 34.1718, -118.1224, -118.02524] #small box I 210 arcadia 

    if delaytype == 'Polynomial': theta = matrix([0.0, 0.0, 0.0, 0.15])
    g = osm_network(theta, delaytype, box)  
    g.visualize(True)
    d.draw(g, nodes = False)
    n = g.numlinks
    print n
    l, x = ue.solver(g, update=True, full=True)
    t1=time.time()
    print 'time for computation for '+str(t1-t0)
    d.draw_ODs(g, 205)
    d.draw_delays(g)

    print max(mul(l,g.get_slopes()))
    print 'cost UE:', sum([link.delay*link.flow for link in g.links.values()])
    #l2, x2 = ue.solver(g, update=True, full=True, SO=True)
    #print (max(mul(l2,g.get_slopes())))
    #print ('cost SO:', sum([link.delay*link.flow for link in g.links.values()]))

def main():
    test2('Polynomial')
    #test_osm_box('Polynomial')
    

if __name__ == '__main__':
    main()