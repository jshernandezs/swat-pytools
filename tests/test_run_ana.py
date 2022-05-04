#!/usr/bin/env python3
import os
import time
import subprocess as sp
from swat_utilities.swat_config import ModelSetup

model_file_path = os.path.abspath('../resources/Models/Test.zip')
swat_model = ModelSetup(model_file_path)

swat_model.swat_exec_name = 'SWAT_Rev681'
swat_model.swat_dir = os.path.abspath('../resources')
swat_model.verbose_swat = True

swat_model.prepare_swat()

start_time = time.time()
swat_model.run_swat()
print("--- Model execution time: {:.5f} seconds ---".format((time.time() - start_time)))

output_dir = '../resources/swat_output/Test_AnaMaria'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

sp.run(['mv', '/tmp/output_swat/New_SWAT/output.rch', os.path.abspath(output_dir)], check=True)
sp.run(['mv', '/tmp/output_swat/New_SWAT/output.sub', os.path.abspath(output_dir)], check=True)
swat_model.remove_swat()
