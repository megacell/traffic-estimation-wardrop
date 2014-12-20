'''
Created on Apr 23, 2014

@author: jeromethai
'''
import scipy.linalg as sla
import scipy.io as sio
import numpy as np
from cvxopt import matrix
import numpy.random as ra
import networkx as nx
import math


def place_zeros(M, tol=1e-13):
    """Replace entries in M less than tol by 0.0"""
    for i in range(M.size[0]):
        for j in range(M.size[1]):
            if M[i,j] < tol: M[i,j] = 0.0
    return M


def find_basis(M):
    """Find the indices of the columns of M that form a basis or range(M)"""
    p,l,u = sla.lu(M)
    ind = [i for i in range(u.shape[0]) if u[i,i] != 0.0]
    if u[i,i] == 0:
        for j in range(i+1,u.shape[1]):
            if u[i,j] != 0.0: ind.append(j); break
    return ind


def bisection(F, f, left, right, tol=1e-8):
    """Use bisection to find x such that F(x)=f
    we suppose F is strictly increasing"""
    l, r = left, right
    while r-l>tol:
        if F((l+r)/2.0) < f:
            l = (l+r)/2.0
        else:
            r = (l+r)/2.0
    return (l+r)/2.0


def save_mat(Ms, names, filename):
    """Save matrices in matlab format
    
    Parameters
    __________
    Ms: list of matrices
    names: list of names of the matrices in matlab format
    """
    dict = {names[k] : np.array(matrix(Ms[k])) for k in range(len(Ms))}
    sio.savemat(filename + '.mat', mdict=dict)


def add_noise(A, a, tol=0.1):
    """add gaussian noise to entries of A that are > tol"""
    m,n = A.size
    M = matrix(0.0, (m,n))
    for i in range(m):
        for j in range(n):
            M[i,j] = A[i,j]
            if M[i,j] > tol: M[i,j] = ra.normal(A[i,j], a*A[i,j])
    return M


def create_networkx_graph(graph):
    """Create networkx.DiGraph graph from graph object"""
    G=nx.DiGraph(indedges={})
    G.add_nodes_from(range(1,graph.numnodes))
    G.add_edges_from([(key[0],key[1]) for key in graph.links.keys()])
    i = 0
    for edge in G.edges(): G.graph['indedges'][edge] = i; i+=1
    return G


def sample_box(N, box):
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
    

def sample_line(N, line, scale):
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


def in_box(point, box):
    """Check if point (x,y) is in box (x1,y1,x2,y2)"""
    x1,y1,x2,y2 = box
    return x1<=point[0]<=x2 and y1<=point[1]<=y2


def distance_on_unit_sphere(lat1, long1, lat2, long2):

    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
        
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
        
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
        
    # Compute spherical distance from spherical coordinates.
        
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta, phi)
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
    
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + 
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )

    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    #return 3960.*arc to get in miles
    #return 6373.*arc to get in km
    return arc

def closest_node(lat, lng, nodes):
    Nodes_candidates=[]
    Distances=[]
    latmax=lat+0.2
    latmin=lat-0.2
    lngmax=lng+0.2
    lngmin=lng-0.2
    for i in range(len(nodes)):
        node = nodes[i]
        if is_in_box(latmin, latmax, lngmin, lngmax, node):
            Nodes_candidates.append([node[0], node[1], i+1]) #i+1 because the order is shifted because the first node id is 1 
    for node in Nodes_candidates:
        Distances.append(distance_on_unit_sphere(lat, lng, node[1], node[0]))
    return Nodes_candidates[np.argmin(Distances)][2]

def is_in_box(latmin, latmax, lngmin, lngmax, node):
    lng = node[0]
    lat = node[1]
    if (lat<latmax and lat>latmin and lng<lngmax and lng>lngmin):
        return True
    else:return False

if __name__ == '__main__':
    pass