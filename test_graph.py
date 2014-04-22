'''
Created on Apr 18, 2014

@author: jeromethai
'''

import graph as g

def smallGrid():
    
    grid = g.grid(2, 3, outDown=[1,1,1], outDdelay=[[(1.0, 4.0)], [(1.0, 4.0)], [(1.0, 1.0)]], 
                  inRight=[1,1,0,1,1], inRdelay=[[(3.0, 1.0)], [(2.0, 1.0)], [(0.0, 0.0)], [(2.0, 1.0)], [(3.0, 1.0)]])

    grid.addOD(3, 6, 3.0)
    grid.addOD(3, 4, 3.0)
    grid.addOD(3, 2, 1.0)
    grid.addOD(2, 4, 1.0)
    
    #grid.addOD(10,3) # error
    #grid.addOD(0,3) # error
    #grid.addOD(2,4) # error
    
    grid.addPath([(3,6,1)])
    grid.addPath([(3,2,1), (2,1,1), (1,4,1)])
    grid.addPath([(3,2,1), (2,5,1), (5,4,1)])
    
    #grid.addPath([(3,2,1), (2,5,1)]) # error
    #grid.addPath([(3,2,1), (2,5,1), (5,4,1)]) # error
    
    grid.addPath([(3,6,1), (6,5,1), (5,4,1)])
    grid.addPath([(3,2,1)])
    grid.addPath([(2,1,1), (1,4,1)])
    grid.addPath([(2,5,1), (5,4,1)])
    return grid

if __name__ == '__main__':
    
    print 'Evening rush example'
    print 

    grid = smallGrid()
    g.visualize(grid, True, True, True, True, True)
    
