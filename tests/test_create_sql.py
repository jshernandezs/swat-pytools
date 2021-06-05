#!/usr/bin/env python3
import time

from swat_utilities.sql_utilities import output_to_sqlite

filenames = ['../resources/swat_output/Honeyoey_Model/output.sub',
             '../resources/swat_output/Honeyoey_Model/output.rch']
output_db = '../resources/swat_output/Honeyoey_Model/output_sqlite.db'

start_time = time.time()
output_to_sqlite(filenames, output_db)
print("--- sqlite database creation: {:.5f} seconds ---".format((time.time() - start_time)))
