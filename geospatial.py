'''
Created on Dec 4, 2014

@author: jeromethai
'''
import numpy as np
import psycopg2
import utm

latmax, latmin, longmax, longmin = 34.288621 , 33.9185, -117.646497, -118.394048

def isInBox(lat, long):
    return (lat < latmax and lat > latmin and long < longmax and long > longmin)
    
    
def selectTAZinBox():
    """Select TAZ ids that are in LA county and in the bounding box specified by isInBox
    """
    tazInBox = {}
    conn = psycopg2.connect(dbname='ca-census', user='postgres')
    cur = conn.cursor()
    cur.execute('SELECT ST_AsText(ST_Centroid(geom)), gid, ST_AsText(geom) FROM taz_data WHERE county = \'037\'')
    for row in cur.fetchall():
        tmp = map(float, row[0][6:-1].split())
        lat, long = utm.to_latlon(tmp[0], tmp[1] , 11, 'N') #shapefiles are in the UTM "11N" coordinates
        if isInBox(lat, long): tazInBox[row[1]] = [lat, long, row[2]]# if centroid is in box
    return tazInBox
        
    
def countNodeInTAZ(tazInBox):
    """count number of nodes in each TAZ
    """
    conn = psycopg2.connect(dbname='ca-census', user='postgres')
    cur = conn.cursor()
    for k,v in tazInBox.items():
        cur.execute("SELECT count(*) FROM la_nodes WHERE ST_Contains( ST_GeomFromText('"+ v[2] +"'), geom)")
        v.append(cur.fetchall()[0][0])
    return tazInBox
    

if __name__ == '__main__':
    tazInfo = countNodeInTAZ(selectTAZinBox())
    for k,v in tazInfo.items(): print k,v[3]
   