#!/usr/bin/env python3
import time
import multiprocessing
import pandas as pd

from swat_utilities.sql_utilities import get_time_series, create_connection

def single_run(db_file, ext, variables, sub):    
    with create_connection(db_file) as db:
        res = get_time_series(db, ext, variables, sub)    
    return res


def run(pool_size):
    
    output_db = '../resources/swat_output/Honeyoey_Model/output_sqlite.db'
    subbasins_filename = '../resources/csv_files/subbasins.csv'    
    
    subbasins = pd.read_csv(subbasins_filename, header=None).iloc[:, 0].to_list()
    var_list = [('rch', ['FLOW_OUTcms']),
                ('sub', ['AREAkm2', 'SYLDt_ha', 'ORGNkg_ha', 'ORGPkg_ha', 'NSURQkg_ha', 'SOLPkg_ha', 'SEDPkg_ha'])]
    
    with multiprocessing.Pool(pool_size) as pool:
        output = []
        for ext, variables in var_list:        
            results = [pool.apply_async(single_run, args=(output_db, ext, variables, x)) for x in subbasins]        
            output.append([r.get() for r in results])          
    
    return output

    
if __name__ == '__main__':
    processes = 8 
    multiprocessing.freeze_support()
    start_time = time.time()
    results = run(processes)
    print("--- sqlite database query: {:.5f} seconds ---".format((time.time() - start_time)))
