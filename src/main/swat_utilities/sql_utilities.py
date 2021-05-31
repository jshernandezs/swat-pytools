#!/usr/bin/env python3
"""
Functions for interacting with SWAT outputs and sqlite
"""
import csv
import os
import sqlite3
import math
import pandas as pd
import datetime as dt
from sqlite3 import Error


def output_to_sqlite(filename_list, output_db, bufsize=100000):
    conn = create_connection(output_db)    
    for file in filename_list:
        ext = file[-3:]
        create_database(conn, file, ext, bufsize)    
    conn.close()

        
def get_time_series(db, ext, variables, unit):
    
    if ext == 'rch':
        table = 'output_rch'
        column = 'RCH'
        time = 'MO, DA, YR'
    elif ext == 'sub':
        table = 'output_sub'
        column = 'SUB'
        time = 'MO, DA, YR'
        variables = ['AREAkm2'] + variables
    else:
        print('Extension not available yet')
        return
    
    cur = db.cursor()
    sql_s = 'SELECT {:s},{:s} FROM {:s} WHERE {:s}=?'.format(time, ','.join(variables), table, column)
    cur.execute(sql_s, (unit,))    
    rows = cur.fetchall()
    
    time = [dt.datetime(year=row[2], month=row[0], day=row[1]) for row in rows]
    values = [row[3:] for row in rows]
    
    d = {variable: [x[i] for x in values] for i, variable in enumerate(variables)}   
    ts = pd.DataFrame(data=d, index=time)
  
    return ts
    

def create_database(db, filename, extension, max_lines=100000):
    if extension  == 'rch':
        table = 'output_rch'
        max_col1 = 44
        max_col2 = 44
        ns = 12
    elif extension == 'sub':
        table = 'output_sub'
        max_col1 = 40
        max_col2 = 41
        ns = 10
    else:
        print('Extension not defined yet.')
        return
    
    with open(os.path.abspath(filename), 'r') as f:
        
        head = [next(f) for x in range(9)]
        raw_header = head[8]
        part1 = raw_header[:max_col1].split()
        part2 = raw_header[max_col1:]
        
        nvar = int(math.ceil((len(part2)-1)/ns))
        part2 = [part2[int((i-1)*ns):int(i*ns)].strip() for i in range(1,nvar+1)]            
        part2 = ['_'.join(x.split('/')) for x in part2]
        part2 = [''.join(x.split('#')) for x in part2]
        part2 = ['_'.join(x.split()) for x in part2]
        
        header = part1 + part2
        
        # generate CREATE TABLE statement
        sql_s = ['CREATE TABLE IF NOT EXISTS {:s} ('.format(table),
                 '{:s} integer,'.format(part1[0]),
                 '{:s} text,'.format(part1[1]),
                 '{:s} integer,'.format(part1[2]),
                 '{:s} integer,'.format(part1[3]),
                 '{:s} integer,'.format(part1[4]),
                 '{:s} real,'.format(part1[5])]
        
        for i, item in enumerate(part2):
            
            if i < (len(part2) - 1):
                sql_s.append('{:s} real,'.format(item))
            else:
                sql_s.append('{:s} real'.format(item))
        
        sql_s.append(');')
        
        # create table
        if db is not None:
            create_table(db, ''.join(sql_s))
        else:
            print("Error! cannot create the database connection.")
            return
        
        reader = csv.reader(f)
        
        c = 1
        buf = []
        
        for line in reader:                
            
            if extension == 'rch':
                part1 = [line[0][5:10], line[0][10:19], line[0][19:23], line[0][23:26], line[0][26:31], line[0][31:41]]                
            elif extension == 'sub':
                part1 = [line[0][6:11], line[0][11:19], line[0][19:23], line[0][23:26], line[0][26:31], line[0][31:41]]
                part1[-1] = '0' + part1[-1]
            
            part1 = [x.strip() for x in part1]               
            part2 = line[0][max_col2:]
            part2 = [part2[int((i-1)*ns):int(i*ns)].strip() for i in range(1,nvar+1)]        
            line_list = part1 + part2
            
            buf.append(tuple(line_list))
            c += 1
            
            if (c > max_lines):                    
                create_records(db, table, header, buf)
                c = 1
                buf = []
        
        if len(buf) > 0:
            create_records(db, table, header, buf)
                    

# definitions for linking sqlite and python; source: https://www.sqlitetutorial.net/sqlite-python/
def create_connection(db_file, time=5):
    """ create a database connection to the SQLite database
        specified by db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file, timeout=time)
    except Error as e:
        print(e)
   
    return conn


def create_table(conn, create_table_sql):
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)
   
     
def create_records(conn, table, columns, list_values):
    
    tokens = ','.join(['?']*len(list_values[0]))
    
    sql = 'INSERT INTO {:s}({:s}) VALUES ({:s})'.format(table, ','.join(columns), tokens)
    cur = conn.cursor()
    cur.executemany(sql, list_values)
    conn.commit()
    
    return cur.lastrowid
