# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 16:31:03 2014

@author: hugo
"""
import numpy as np
from util import distance_on_unit_sphere


class Manips:
    def __init__(self):
        self.List_TAZ = np.genfromtxt('LA_medium_data/Description_TAZ.csv', delimiter=',', skip_header = 1)
        self.dict_TAZ ={}
        for i in range(len(self.List_TAZ)):
            self.dict_TAZ[self.List_TAZ[i][0]] = i 
        	
        self.List_TAZ_ids = self.List_TAZ[:,0].astype(int)
              
    def Read_TAZ_from_csv(self):
        flows = []
        i=0
        with open('LA_medium_data/CTPP_LA.csv', 'rb') as f:
            for line in f:
			l = line.split('\",\"')
			i=i+1
			tazi = int( l[0].split(',')[0].split(' ')[1]) # parsed TAZ_id
			tazi = int(-20000+tazi/1000)
			tazj = int(l[1].split(',')[0].split(' ')[1]) # parsed TAZ_id
			tazj = int(-20000+tazj/1000)	
			if (tazi in self.List_TAZ_ids and tazj in self.List_TAZ_ids):
				f = float(l[4].replace(',','')) + float(l[6].replace(',',''))/2.0 + float(l[8].replace(',',''))/3.0 + float(l[10].replace(',',''))/4.0 + float(l[12].replace(',',''))/5.0 + float(l[14].replace(',',''))/6.0     # total vehicle-trips, we divide when people are a lot in cars
				if (f!=0 and tazi !=tazj):				
					flows.append([tazi, tazj, float(f)])
        return flows
    
    def Dist_between_TAZ(self, tazi, tazj):
        return distance_on_unit_sphere(self.List_TAZ[self.dict_TAZ[tazi]][1], self.List_TAZ[self.dict_TAZ[tazi]][2], self.List_TAZ[self.dict_TAZ[tazj]][1], self.List_TAZ[self.dict_TAZ[tazj]][2])

def dist_to_I210box(lat, lng, box = [34.124918, 34.1718, -118.1224, -118.02524]):
    if lat > box[1]:
        if lng < box[2]: return distance_on_unit_sphere(lat, lng, box[1], box[2])
        if lng > box[3]: return distance_on_unit_sphere(lat, lng, box[1], box[3])
        else : return distance_on_unit_sphere(lat, lng, box[1], lng)
    if lat < box[0]:
        if lng < box[2]: return distance_on_unit_sphere(lat, lng, box[0], box[2])
        if lng > box[3]: return distance_on_unit_sphere(lat, lng, box[0], box[3])
        else : return distance_on_unit_sphere(lat, lng, box[0], lng)
    elif lng < box[2]: return distance_on_unit_sphere(lat, lng, lat, box[2])
    elif lng > box[3]: return distance_on_unit_sphere(lat, lng, lat, box[3])
    else : return 0

def Is_in_I210box(lat, lng, type = 'box'):
    if type == 'box' :box = [34.124918, 34.1718, -118.1224, -118.02524]
    if type == 'medium' : box = [34.081133, 34.237951, -118.249853, -117.893484]
    return (lat < box[1] and lat > box[0] and lng < box[3] and lng > box[2])

'''
m=Manips()
ODs = np.asarray(m.Read_TAZ_from_csv())
'''