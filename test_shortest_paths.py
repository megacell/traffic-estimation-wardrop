'''
Created on Aug 7, 2014

@author: jeromethai
'''

import shortest_paths as sh
from generate_graph import los_angeles


def test1():
    g = los_angeles()[0]
    dist, next = sh.Dijkstra(g, 1, [17])
    #print sh.get_path(17, 1, next)
    sink = 20
    sources = [od[0] for od in g.ODs.keys() if od[1]==sink]
    #print sources
    sources = [6]
    As = sh.mainKSP(g, sources, sink, 10)
    for s in sources:
        print len(As[s]), As[s]
        for i in range(len(As[s])):
            for j in range(i+1,len(As[s])):
                if As[s][i]==As[s][j]: print 'same list'
            
    

def main():
    test1()


if __name__ == '__main__':
    main()