#!/usr/bin/env python3
"""
Write SWAT txt files based on user-defined crop rotations and operations
"""

import os
import platform
import time
import pandas as pd

import swat_utilities.rotation_factory as builder
from swat_utilities.swat_config import ModelSetup

# SETTINGS

# 1) Scenario settings

opt1 = 2  # Select rotation
opt2 = 1  # Select system
opt3 = 3  # Select application rate
opt4 = 2  # Select cover crop
opt5 = 1  # Select tile drainage option
opt6 = 1  # Select filter strip option

start_fs = {'MONTH': 1, 'DAY': 1, 'IYEAR': 1998}
subbasins_filename = '../resources/csv_files/subbasins.csv'

# 2) Model settings

model_file_path = os.path.abspath('../resources/Models/SWAT_Model.zip')
if platform.system() == 'Linux':
    swat_version = 'SWAT_Rev670'
elif platform.system() == 'Darwin':
    swat_version = 'SWAT_Rev622_macOS'
else:
    swat_version = None
    print('Windows version is not supported yet')

# MAIN SCRIPT STARTS HERE

# 1) Generate schedule variables

rotation, params, filter_strip = builder.set_rot(opt1, opt2, opt3, opt4, opt5, opt6, start_fs)
params.update({'NBYR': [12, 'replace', 'cio'],
               'NYSKIP': [2, 'replace', 'cio'],
               'IPRINT': [2, 'replace', 'cio'],
               'ICALEN': [0, 'replace', 'cio']})

# 2) Prepare SWAT model

comb = [opt1, opt2, opt3, opt4, opt5, opt6]
subbasins = pd.read_csv(subbasins_filename, header=None).iloc[:, 0].to_list()

swat_model = ModelSetup(model_file_path)
swat_model.swat_exec_name = swat_version
swat_model.verbose_swat = True
swat_model.param = [params]
swat_model.drainmod = True
swat_model.mgt_table = [rotation]
swat_model.ops_table = [filter_strip]
swat_model.subbasins = [subbasins]
swat_model.subbasins_mgt = [subbasins]
swat_model.subbasins_ops = [subbasins]
swat_model.new_model_name = '{:d}{:02d}{:d}{:d}{:d}{:d}'.format(*comb)

start_time = time.time()
swat_model.prepare_swat()
print("--- Pre-processing time was {:.2f} minutes ---".format((time.time() - start_time) / 60))

# 3) Run SWAT model
start_time = time.time()
swat_model.run_swat()
print("--- Simulation time was {:.2f} minutes ---".format((time.time() - start_time) / 60))
