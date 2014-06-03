'''
Created on Apr 18, 2014

@author: jeromethai
'''

import Graph as g


def small_example():
    
    graph = g.Graph('Small example graph')
    
    graph.add_node((0,1))
    graph.add_node((0,-1))
    graph.add_node((2,0))
    graph.add_node((4,0))
    graph.add_node((6,0))
    
    graph.add_link(1, 3, 1, delayfunc=g.create_delayfunc('Affine',(1.0, 0.0)))
    graph.add_link(2, 3, 1, delayfunc=g.create_delayfunc('Affine',(2.0, 0.0)))
    graph.add_link(3, 4, 1, delayfunc=g.create_delayfunc('Affine',(2.0, 1.0)))
    graph.add_link(3, 4, 2, delayfunc=g.create_delayfunc('Affine',(1.0, 2.0)))
    graph.add_link(4, 5, 1, delayfunc=g.create_delayfunc('Affine',(1.0, 0.0)))
    
    graph.add_od(1, 5, 2.0)
    graph.add_od(2, 5, 3.0)
    
    graph.add_path([(1,3,1), (3,4,1), (4,5,1)])
    graph.add_path([(1,3,1), (3,4,2), (4,5,1)])
    graph.add_path([(2,3,1), (3,4,1), (4,5,1)])
    graph.add_path([(2,3,1), (3,4,2), (4,5,1)])
    
    return graph


def small_grid():
    
    grid = g.create_grid(2, 3, outdown=[1,1,1], outddelay=[[(1.0, 4.0)], [(1.0, 4.0)], [(1.0, 1.0)]], 
                  inright=[1,1,0,1,1], inrdelay=[[(3.0, 1.0)], [(2.0, 1.0)], [(0.0, 0.0)], [(2.0, 1.0)], [(3.0, 1.0)]], delaytype='Affine')
    
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


def small_grid2():
    grid = g.create_grid(2, 3, outdown=[1,1,1], outddelay=[[(1.0, 4.0)], [(1.0, 4.0)], [(1.0, 1.0)]], 
                  inright=[1,1,0,1,1], inrdelay=[[(3.0, 1.0)], [(2.0, 1.0)], [(0.0, 0.0)], [(2.0, 1.0)], [(3.0, 1.0)]], delaytype='Affine')
    
    grid.add_link(5, 2, 1, delayfunc=g.AffineDelay(1.0, 4.0))
    
    grid.add_od(3, 6, 1.0)
    grid.add_od(3, 4, 1.0)
    grid.add_od(3, 2, 1.0)
    grid.add_od(2, 4, 4.0)
    
    grid.add_path([(3,6,1)])
    grid.add_path([(3,2,1), (2,1,1), (1,4,1)])
    grid.add_path([(3,2,1), (2,5,1), (5,4,1)])
    grid.add_path([(3,6,1), (6,5,1), (5,4,1)])
    grid.add_path([(3,6,1), (6,5,1), (5,2,1), (2,1,1), (1,4,1)])
    grid.add_path([(3,2,1)])
    grid.add_path([(3,6,1), (6,5,1), (5,2,1)])
    grid.add_path([(2,1,1), (1,4,1)])
    grid.add_path([(2,5,1), (5,4,1)])
    return grid


def main():
    #graph = small_example()
    graph = small_grid()
    graph.visualize(True, True, True, True, True)    


if __name__ == '__main__':
    main()
    
    
