'''
Created on Apr 18, 2014

@author: jeromethai
'''

import Graph as g

def small_grid():
    
    grid = g.create_grid(2, 3, outdown=[1,1,1], outddelay=[[(1.0, 4.0)], [(1.0, 4.0)], [(1.0, 1.0)]], 
                  inright=[1,1,0,1,1], inrdelay=[[(3.0, 1.0)], [(2.0, 1.0)], [(0.0, 0.0)], [(2.0, 1.0)], [(3.0, 1.0)]])
    
    grid.add_link(5, 2, 1, delayfunc=g.AffineDelay(1.0, 4.0))
    
    grid.add_od(3, 6, 3.0)
    grid.add_od(3, 4, 3.0)
    grid.add_od(3, 2, 1.0)
    grid.add_od(2, 4, 1.0)
    
    #grid.add_od(10,3) # error
    #grid.add_od(0,3) # error
    #grid.add_od(2,4) # error
    
    # paths for OD (3,6)
    grid.add_path([(3,6,1)])
    
    #paths for OD (3,4)
    grid.add_path([(3,2,1), (2,1,1), (1,4,1)])
    grid.add_path([(3,2,1), (2,5,1), (5,4,1)])
    grid.add_path([(3,6,1), (6,5,1), (5,4,1)])
    grid.add_path([(3,6,1), (6,5,1), (5,2,1), (2,1,1), (1,4,1)])
    
    #grid.add_path([(3,2,1), (2,5,1)]) # error
    #grid.add_path([(3,2,1), (2,5,1), (5,4,1)]) # error
    
    # paths for OD (3,2)
    grid.add_path([(3,2,1)])
    grid.add_path([(3,6,1), (6,5,1), (5,2,1)])
    
    # paths for OD (2,4)
    grid.add_path([(2,1,1), (1,4,1)])
    grid.add_path([(2,5,1), (5,4,1)])
    return grid


def main():
    print 'Evening rush example'; print 
    grid = small_grid()
    g.visualize(grid, True, True, True, True, True)    


if __name__ == '__main__':
    main()
    
    
