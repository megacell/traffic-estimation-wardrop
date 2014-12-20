import graph as g
import numpy as np
from cvxopt import matrix
import scipy.io as sio

class grid:
        
    def __init__(self, nodes_on_a_side, nodes_on_other_side, taz_size):
        self.nodes_on_a_side = nodes_on_a_side
        self.nodes_on_other_side = nodes_on_other_side
        self.taz_size = taz_size
        self.nodes_list, self.links_list = self.make_grid(self.nodes_on_a_side, self.nodes_on_other_side)
        print self.nodes_list        
        #self.ODs = self.create_TAZ_ODs_for_grid()
        self.ODs=[[self.nodes_on_a_side*self.nodes_on_other_side, self.nodes_on_a_side*self.nodes_on_other_side - 1, 5.0]]            
        print self.ODs        
        self.graph = g.create_graph_from_list(self.nodes_list, self.links_list, 'Polynomial', self.ODs, 'Map of L.A.')
    
    def create_TAZ_ODs_for_grid(self):
        ODsTaz=[]
        data = sio.loadmat('los_angeles_data_2.mat')   
        ODs = data['ODs4']
        for i in range(1):
        #for i in range(len(ODs)):
            startnode, endnode, od_value = ODs[i][0], ODs[i][1], ODs[i][2]
            starttaz = self.Which_Taz(startnode)
            endtaz = self.Which_Taz(endnode)
            print self.taz_size
            print startnode
            print starttaz
            ODsTaz.append([(starttaz[0]*self.taz_size + 1)*self.nodes_on_a_side +1 + starttaz[1]*self.taz_size, (endtaz[0]*self.taz_size + 1)*self.nodes_on_a_side+1+endtaz[1]*self.taz_size, od_value/4])
            ODsTaz.append([((starttaz[0]+1)*self.taz_size - 1)*self.nodes_on_a_side+1 + starttaz[1]*self.taz_size, ((endtaz[0]+1)*self.taz_size - 1)*self.nodes_on_a_side+1+endtaz[1]*self.taz_size, od_value/4])
            ODsTaz.append([(starttaz[0]*self.taz_size +1)*self.nodes_on_a_side + self.taz_size-1 + starttaz[1]*self.taz_size, (endtaz[0]*self.taz_size +1)*self.nodes_on_a_side + self.taz_size-1+endtaz[1]*self.taz_size, od_value/4])
            ODsTaz.append([((starttaz[0]+1)*self.taz_size - 1)*self.nodes_on_a_side + self.taz_size-1 + starttaz[1]*self.taz_size, ((endtaz[0]+1)*self.taz_size - 1)*self.nodes_on_a_side + self.taz_size-1+endtaz[1]*self.taz_size, od_value/4])
        return ODsTaz
        
    def Which_Taz(self, node):
        T=[[2,(2,3)],[8,(4,3)],[9,(6,3)],[10,(7,3)],[4,(0,2)],[5,(1,2)],[6,(2,2)],[7,(2,2)],[17,(4,2)],[18,(6,2)],[19,(7,2)],[11,(0,1)],[22,(2,1)],[24,(4,1)],[26,(6,1)],[27,(7,1)],[20,(0,0)],[29,(1,0)],[30,(3,0)],[31,(3,0)],[25,(6,0)],[33,(6,0)]]
        for t in T:
            if t[0] == node: return t[1]
        
    def assign_ODs_to_TAZ(self, OD_value, TAZ_id_i, TAZ_id_j):
        
        return 1
    
    def make_grid(self, nodes_on_a_side, nodes_on_other_side):
        length_of_a_link = 1.0
        speed_limit=32
        nodes_list = []
        for i in range(nodes_on_a_side):
            for j in range(nodes_on_other_side):
                if (i+j==0): nodes_list.append([i*length_of_a_link, j*length_of_a_link, nodes_on_a_side*nodes_on_other_side])
                else: nodes_list.append([i*length_of_a_link, j*length_of_a_link, i*nodes_on_a_side + j])
        
        links_list=[]
        for i in range(nodes_on_a_side):
            for j in range(nodes_on_other_side):
                if (i!=nodes_on_a_side-1 and i+j!=0):
                    links_list.append([i*nodes_on_a_side + j, (i+1)*nodes_on_a_side + j, 1, length_of_a_link/speed_limit, 1/5.0])
                if (i!=0 and i*nodes_on_a_side + j!=1 and i*nodes_on_a_side + j!=nodes_on_a_side):
                    links_list.append([i*nodes_on_a_side + j, (i-1)*nodes_on_a_side + j, 1, length_of_a_link/speed_limit, 1/5.0])
                if (j!=nodes_on_other_side-1 and i+j!=0):
                    links_list.append([i*nodes_on_a_side + j, i*nodes_on_a_side +j+1, 1, length_of_a_link/speed_limit, 1/5.0])
                if (j!=0 and i*nodes_on_a_side + j!=1 and i*nodes_on_a_side + j!=nodes_on_a_side):
                    links_list.append([i*nodes_on_a_side + j, i*nodes_on_a_side + j-1, 1, length_of_a_link/speed_limit, 1/5.0])
                if (i*nodes_on_a_side +j==1):
                    links_list.append([1, nodes_on_a_side*nodes_on_other_side, 1, length_of_a_link/speed_limit, 1/5.0])
                if (i*nodes_on_a_side +j==nodes_on_a_side):
                    links_list.append([nodes_on_a_side, nodes_on_a_side*nodes_on_other_side, 1, length_of_a_link/speed_limit, 1/5.0])
                if (i+j==0):
                    links_list.append([nodes_on_a_side*nodes_on_other_side, 1, 1, length_of_a_link/speed_limit, 1/5.0])
                    links_list.append([nodes_on_a_side*nodes_on_other_side, nodes_on_a_side, 1, length_of_a_link/speed_limit, 1/5.0])
    
        tmp = links_list
        links_list=[]
        theta = matrix([0,0.15])
        degree = len(theta)
        for startnode, endnode, route, ffdelay, slope in tmp:
            coef = [ffdelay*a*b for a,b in zip(theta, np.power(slope, range(1,degree+1)))]
            links_list.append((startnode, endnode, route, ffdelay, (ffdelay, slope, coef))) 
            
        return nodes_list, links_list
'''
        ODs=[[nodes_on_a_side**2, nodes_on_a_side**2-1, 5.0]]
        print ODs
        return g1
'''