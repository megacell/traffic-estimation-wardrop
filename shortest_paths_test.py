'''
Created on Aug 7, 2014

@author: jeromethai
'''

import shortest_paths as sh
from generate_graph import los_angeles


def test1():
    g = los_angeles()[0]
    print g.links.values()[23].delay
    sink = 20
    sources = [od[0] for od in g.ODs.keys() if od[1]==sink]
    As = sh.mainKSP(g, sources, sink, 10)
    for s in sources:
        print len(As[s]), As[s]
    

def main():
    test1()


if __name__ == '__main__':
    main()