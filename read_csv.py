# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 11:49:41 2014

@author: hugo
"""
import numpy as np
from util import is_in_box

def create_nodes_links_from_box(latmin=33.0, latmax=36.0, lngmin=-120.0, lngmax=-116.0):
    data = np.genfromtxt('Data/Network/CSV/LA_small_box_detailed/nodes_la_exported2.csv', delimiter = ';', skiprows=1)
    Nodes_table=[]
    osmId_2_gId = {}
    k = 1
    for i in range(len(data)):
        node=[data[i][0], data[i][1], data[i][2]]
        if is_in_box(latmin, latmax, lngmin, lngmax, node):
            Nodes_table.append([node[0], node[1]])    
            osmId_2_gId[int(node[2])] = k # osmId_2_gId[int(osmId)] = gId
            k += 1
    print str(osmId_2_gId[1687115183])+ ' this is the dest'
    data = np.genfromtxt('Data/Network/CSV/LA_small_box_detailed/links_la_exported.csv', delimiter = ',', skiprows=1)
    link_table=[]
    for i in range(len(data)):
        sn, en=int(data[i][1]), int(data[i][2])
        if (sn in osmId_2_gId.keys() and en in osmId_2_gId.keys()):
            startnode = osmId_2_gId[sn]
            endnode = osmId_2_gId[en]
            length = data[i][3]
            cap = data[i][5]
            freespeed = data[i][6]
            ff_d=length/freespeed
            #slope= 2000 / cap
            slope= 1 / cap
            u=[startnode,endnode,1,ff_d,slope]
            link_table.append(u)
    return Nodes_table, link_table
         
def kill_double_links(links_list):
    links_without_double=[]
    counter=0
    for i in range(len(links_list)):
        link=links_list[i]
        startnode,endnode=link[0], link[1]
        flag=True
        for j in range(len(links_list)-i-1):
            link_j=links_list[j+i+1]
            if (startnode==link_j[0] and endnode == link_j[1]): flag=False
        if (flag==True): links_without_double.append(link)
        else: counter+=1
    print str(counter)+' links removed'
    return links_without_double