#!/usr/bin/env python3
"""
Run a single SWAT model

"""
import os
import platform
import subprocess as sp
import numpy as np
import pandas as pd
from swat_utilities.swat_config import ModelSetup
from swat_utilities.performance import Metric
from swat_utilities.performance import read_obs_file
from swat_utilities.hit import HydrologicIndex as Hi

# Step 1: SWAT model settings
if platform.system() == 'Linux':
    swat_version = 'SWAT_Rev622'
elif platform.system() == 'Darwin':
    swat_version = 'SWAT_Rev622_macOS'
else:
    swat_version = None
    print('Windows version is not supported yet')

model_file_path = os.path.abspath('../resources/Models/Honeyoy_Model.zip')
swat_model = ModelSetup(model_file_path)
swat_model.swat_exec_name = swat_version
swat_model.verbose_swat = True

subbasins_filename = '../resources/csv_files/subbasins.csv'
subbasins = pd.read_csv(subbasins_filename, header=None).iloc[:, 0].to_list()

# Step 2: Parameters to change

params = {'ICALEN': [1, 'replace', 'cio'],
          'BIOMIX': [0.22, 'replace', 'mgt'],
          'CN2': [-0.21, 'multiply', 'mgt'],
          'CANMX': [1.67, 'replace', 'hru'],
          'ESCO': [0.70, 'replace', 'hru'],
          'EPCO': [0.0059, 'replace', 'hru'],
          'GW_DELAY': [6.11, 'replace', 'gw'],
          'ALPHA_BF': [0.83, 'replace', 'gw'],
          'GWQMN': [438, 'replace', 'gw'],
          'GW_REVAP': [0.16, 'replace', 'gw'],
          'REVAPMN': [438, 'replace', 'gw'],
          'RCHRG_DP': [0.50, 'replace', 'gw'],
          'CH_N2': [0.12, 'replace', 'rte'],
          'CH_K2': [6.45, 'replace', 'rte'],
          'SOL_AWC': [-0.21, 'multiply', 'sol'],
          'SURLAG': [1.10, 'replace', 'bsn']}

swat_model.param = [params]
swat_model.subbasins = [subbasins]

# Step 3: Prepare text input files and run SWAT
swat_model.prepare_swat()
swat_model.run_swat()

# Step 4: Post-processing
swat_model.read_output()
swat_model.simulated('watout.dat', plot=False)
summary_model = swat_model.__dict__

# Step 5: Compute performance metrics and make plot
observed_file = '../resources/Observed/Honeyoy_cal.csv'
simulation = Metric(observed_file, swat_model.sim_ts['FLOWm^3/s'])
simulation.get_all()
simulation.plot()
results = simulation.__dict__

# Step 6: Compute hydrologic indicators
observed = read_obs_file(observed_file)
simulated = swat_model.sim_ts['FLOWm^3/s']
obs_hit = Hi(observed)
sim_hit = Hi(simulated, obs_file=observed_file)

index_obs_dict = obs_hit.get_all_hit_indices()
index_sim_dict = sim_hit.get_all_hit_indices()

error = {}
for key in list(index_obs_dict.keys()):
    obs = np.array(index_obs_dict[key])
    sim = np.array(index_sim_dict[key])
    error[key] = abs(obs - sim)/obs * 100

# Step 7: Remove model folder
output_dir = '../resources/swat_output/Honeyoey_Model'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

sp.run(['mv', '/tmp/output_swat/New_SWAT/output.rch', os.path.abspath(output_dir)], check=True)
sp.run(['mv', '/tmp/output_swat/New_SWAT/output.sub', os.path.abspath(output_dir)], check=True)
swat_model.remove_swat()
