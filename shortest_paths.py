'''
Created on Aug 7, 2014

@author: jeromethai
'''


import numpy as np


def Dijkstra(graph, sink, sources=None):
    """Find the shortest path in ffdelays to sink from every other vertex
    Stops when the shortest path form sources to sink have been found
    (see http://en.wikipedia.org/wiki/Dijkstra_algorithm)
    
    Return value:
    -------------
    dist: dist[u] = distance from u to sink
    next: next[u] = next node in the shortest path from u to sink
    """
    n = graph.numnodes
    dist = {k+1:np.inf for k in range(n)}
    next = {k+1:None for k in range(n)}
    dist[sink] = 0.0
    S, Q = list(sources[:]), range(1,n+1)
    while len(Q)>0:
        min = np.inf
        for v in Q:
            if dist[v] < min: u=v; min = dist[v]
        if min == np.inf: return dist, next
        for s in S:
            if u == s: S.remove(u)
        if len(S)==0: return dist,next
        Q.remove(u)
        for link in graph.nodes[u].inlinks.values():
            v, alt = link.startnode, dist[u] + link.ffdelay
            if alt < dist[v]: dist[v] = alt; next[v] = u
    return dist, next


def get_path(source, sink, next):
    """Get shortest path from source to sink from next (given by Dijkstra)
    
    Return value:
    -------------
    Shortest path: list of nodes
    """
    u, path = source, [source]
    while u != sink: u=next[u]; path.append(u)
    return path


def mainKSP(graph, sources, sink, K):
    """Find the K-shortest paths from sources to sink
    
    Return value:
    -------------
    As: dictionary s.t. As[s]=[K-shortest paths from s to sink for s in sources]
    """
    dist, next = Dijkstra(graph, sink, sources)
    A0s = {s:get_path(s, sink, next) for s in sources}
    return {s : YenKSP(graph, s, sink, K, A0s[s]) for s in sources}
            
            
def YenKSP(graph, source, sink, K, A0):
    """"Find the k-shortest paths from source to sink
    A0: initialization with the shortest path from source to sink
    {see http://en.wikipedia.org/wiki/Yen's_algorithm}
    """
    A, B, costs, j, tmp, k2 = [A0], {}, {}, 0, {}, 0
    print 'start YenKSP with A0 =', A
    for k in range(K-1):
        print '################## ITERATION k =', k
        for i in range(len(A[k2])-1):
            spurNode, rootPath = A[k2][i], A[k2][:i+1]
            print '########## spurNode, rootPath =', spurNode, rootPath
            costRootPath = 0
            for l in range(i):
                costRootPath += graph.links[(A[k2][l],A[k2][l+1],1)].ffdelay
            #print 'costRootPath =', costRootPath
            for p in A:
                if rootPath == p[:i+1]:
                    print 'Found same rootpath: ', p
                    if (p[i],p[i+1],1) not in tmp.keys():
                        print 'Remove ', (p[i],p[i+1])
                        link = graph.links[(p[i],p[i+1],1)]
                        tmp[(p[i],p[i+1],1)] = link.ffdelay #save p.edge(i, i + 1)
                        link.ffdelay = np.inf #remove p.edge(i, i + 1)
            for node in rootPath:
                print 'Make node unreachable: ', node 
                for link in graph.nodes[node].inlinks.values():
                    if (link.startnode,node,1) not in tmp.keys():
                        #print 'Remove ', (link.startnode,node)
                        tmp[(link.startnode,node,1)] = link.ffdelay #save edge
                        link.ffdelay = np.inf #remove edge
            #print sink,spurNode
            #print tmp
            dist, next = Dijkstra(graph, sink, [spurNode])
            cost = costRootPath + dist[spurNode]
            if dist[spurNode] < np.inf and cost not in costs.values():
                B[j] = rootPath[:-1] + get_path(spurNode, sink, next)
                costs[j] = cost
                print 'Add path to B with cost: ', B[j], costs[j]
                j += 1
            else:
                print 'Do not path to B', 
            print 'Restore paths and nodes'
            for id,ffdelay in tmp.items(): graph.links[id].ffdelay = ffdelay
            tmp = {}
        if len(B) == 0: break
        min_cost = min(costs.values())
        for key,cost in costs.items():    
            if cost == min_cost:
                print 'Found min cost path of B: ', B[key]
                A.append(B[key]); k2+=1
                print 'Add min cost path of B to A: ', B[key]
                print 'A is now: ', A
                del costs[key]
                del B[key]
                break
        print 'B is now: ', B
    return A
                
            

if __name__ == '__main__':
    pass