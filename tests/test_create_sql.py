#!/usr/bin/env python3
import sys
import os
import time

main_dir = os.path.abspath('../main/')
sys.path.insert(1, main_dir)

from swat_utilities.sql_utilities import output_to_sqlite

filenames = ['../../resources/swat_output/Honeyoy_Model/output.sub',
             '../../resources/swat_output/Honeyoy_Model/output.rch']
output_db = '../../resources/swat_output/Honeyoy_Model/output_sqlite.db'

start_time = time.time()
output_to_sqlite(filenames, output_db)
print("--- sqlite database creation: {:.5f} seconds ---".format((time.time() - start_time)))
