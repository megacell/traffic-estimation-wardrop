'''
Created on Aug 9, 2014

@author: jeromethai
'''

import Waypoints as wp


def test1():
    print wp.box(100, (0.0, 0.0, 1.0, 1.0))
    print wp.line(100, (0.0,0.0,10.0,10.0), 1.0)
    

def main():
    test1()
    

if __name__ == '__main__':
    main()