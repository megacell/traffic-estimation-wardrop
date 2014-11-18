'''
Created on Nov 18, 2014

@author: jeromethai
'''

import psycopg2

def main():
    print 'import sucessful'
    #conn = psycopg2.connect(dbname='template_postgis', user='postgres',
    #    host='localhost', password='jht2115')
    conn = psycopg2.connect(dbname='template_postgis', user='postgres')
    print 'connection successful'
    cur = conn.cursor()
    cur.execute("""SELECT srtext FROM spatial_ref_sys WHERE srid = 36910;""")
    rows = cur.fetchall()
    for row in rows: print "    ", row[1]

if __name__ == '__main__':
    main()