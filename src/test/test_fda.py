#!/usr/bin/env python3
import sys
import os

main_dir = os.path.abspath('../main/')
sys.path.insert(1, main_dir)

from swat_utilities.swat_config import ModelSetup
from swat_utilities.performance import Metric
from swat_utilities.fda_swat import MetricFda

# Step 1: SWAT model settings
modelFilePath = os.path.abspath('../../resources/Models/Honeyoy_Model.zip')
swat_model = ModelSetup(modelFilePath)
swat_model.SWATexecName = 'SWAT_Rev622'

# Step 2: Parameters to change
swat_model.param = {'BIOMIX': [0.10553834, 'replace', 'mgt'],
                    'CN2': [-0.00314644964, 'multiply', 'mgt'],
                    'CANMX': [16.1969501, 'replace', 'hru'],
                    'ESCO': [0.9053373, 'replace', 'hru'],
                    'EPCO': [0.675994726, 'replace', 'hru'],
                    'GW_DELAY': [0.119799926, 'replace', 'gw'],
                    'ALPHA_BF': [0.92920519, 'replace', 'gw'],
                    'GWQMN': [231.704881, 'replace', 'gw'],
                    'GW_REVAP': [0.0733707439, 'replace', 'gw'],
                    'REVAPMN': [739.468178, 'replace', 'gw'],
                    'RCHRG_DP': [0.697223355, 'replace', 'gw'],
                    'CH_N2': [0.0563435942, 'replace', 'rte'],
                    'CH_K2': [97.1313234, 'replace', 'rte'],
                    'SOL_AWC': [-0.0992068353, 'multiply', 'sol'],
                    'SURLAG': [11.8943658, 'replace', 'bsn']}

# Step 3: Prepare text input files and run SWAT
swat_model.prepare_swat()
swat_model.run_swat()

# Step 4: Post-processing
swat_model.read_output()
swat_model.simulated('watout.dat', plot=False)
summary_model = swat_model.__dict__

# Step 5: Compute performance metrics and make plot
observed_file = '../../resources/Observed/Honeyoy.csv'
simulation = Metric(observed_file, swat_model.sim_ts['FLOWm^3/s'])
simulation.get_all()
simulation.plot()
results = simulation.__dict__

# Step 6: FDA metric
obs_file = os.path.abspath('../main/fda_swat/Data/observed.csv')
sim_file = os.path.abspath('../main/fda_swat/Data/simulated.csv')
rds_file = os.path.abspath('../main/fda_swat/Data/observed_fmodel_fast_test.rds')
simulation.write_obs(obs_file)
simulation.write_sim(sim_file)
fda_swat = MetricFda(fda_dir='../main/fda_swat')
fda_swat.fda_fit(obs_file, rds_file, opt_fast=True)
fda_swat.fda_metric(sim_file, rds_file, opt_fast=True)

# Step 7: Remove model folder
swat_model.remove_swat()
