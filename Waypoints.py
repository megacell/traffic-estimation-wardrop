'''
Created on Aug 9, 2014

@author: jeromethai
'''

import numpy.random as ra


class Waypoints:
    """Class Waypoints containing a bbox"""
    def __init__(self, wp={}, N=0, lines={}):
        self.wp = wp # dictionary of all waypoints
        self.N = N # total number of waypoints
        self.lines = lines # total 
        


def box(N, box):
    """Sample uniformly at random N points in a box
    
    Parameters:
    ----------
    N: number of points
    box: (x1,y1,x2,y2) with (x1,y1) lower-left corner, (x2,y2) upper-right corner
    """
    x1,y1,x2,y2 = box
    if x1 > x2: print 'Error: x1 > x2'; return
    if y1 > y2: print 'Error: y1 > y2'; return
    return zip(ra.uniform(x1,x2,N), ra.uniform(y1,y2,N))
    

def line(N, line, scale):
    """Sample uniformly N points on a line with Gaussian noise
    
    Parameters:
    ----------
    N: number of points
    line: (x1,y1,x2,y2) with (x1,y1) one end and (x2,y2) the other end
    scale: standard deviation (spread or width) of the distribution
    """
    x1,y1,x2,y2 = line
    u, s1, s2 = ra.uniform(size=N), ra.normal(scale=scale, size=N), ra.normal(scale=scale, size=N)
    l = [(x1+p*(x2-x1)+t1, y1+p*(y2-y1)+t2) for p,t1,t2 in zip(u,s1,s2)]
    return l



if __name__ == '__main__':
    pass