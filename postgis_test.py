'''
Created on Nov 18, 2014

@author: jeromethai
'''

import psycopg2
import utm # https://pypi.python.org/pypi/utm


def test_utm():
    print utm.from_latlon(51.2, 7.5)

def connect():
    conn = psycopg2.connect(dbname='template_postgis', user='postgres')
    cur = conn.cursor()
    cur.execute("""SELECT srtext FROM spatial_ref_sys WHERE srid = 26910;""")
    #cur.execute("""SELECT srtext FROM spatial_ref_sys;""")
    rows = cur.fetchall()
    for row in rows: print row


def select():
    conn = psycopg2.connect(dbname='ca-census', user='postgres')
    cur = conn.cursor() # creates cursor that establishes connection
    cur.execute("""SELECT * FROM ca_census_data;""")
    rows = cur.fetchall()
    for row in rows: print row


def main():
    select()
    #create_table()
    #connect()
    #test_utm()

if __name__ == '__main__':
    main()