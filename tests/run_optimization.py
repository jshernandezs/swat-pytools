#!/usr/bin/env python3
"""
Generic optimization routine using monthly streamflow and soil water outputs
from SWAT

"""

import os
import time
import pandas as pd
import autograd.numpy as anp
from pymoo.algorithms.moo.unsga3 import UNSGA3
from pymoo.factory import get_sampling, get_crossover, get_mutation, get_reference_directions
from pymoo.optimize import minimize
from swat_utilities.optimization_factory import SWATConfig, SWATProblem, my_callback
#from optimization_factory import SWATConfig, SWATProblem, my_callback
import dask
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
# 


# MODELING SETTINGS ####################################################################################################

output_dir = '../resources/swat_output/Test_AnaMaria'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
X_dir = os.path.abspath('../resources/csv_files/cienaga_int_popu.csv')

# 1) create calibration settings object
swat_model = SWATConfig()

# 2) General settings
swat_model.model_file = os.path.abspath('../resources/Models/model_cienaga.zip')
swat_model.simulation_months = os.path.abspath('../resources/csv_files/cienaga_simulation_months_hru.csv')
swat_model.soils_awc = os.path.abspath('../resources/csv_files/cienaga_awc_soils.csv')
swat_model.hid_treshold = os.path.abspath('../resources/csv_files/rancho_hid_dro_treshold.csv')
swat_model.opt_hrus = os.path.abspath('../resources/csv_files/cienaga_hru_opt.csv')
swat_model.opt_sub = os.path.abspath('../resources/csv_files/cienaga_sub_opt.csv')
swat_model.loc_const_1 = os.path.abspath('../resources/csv_files/cienaga_loc_const.csv') 

swat_model.swat_exec_name = 'SWAT_Rev681'
swat_model.temp_dir = '/tmp/swat_runs'
swat_model.temp_run_dir = '/tmp/output_swat'
swat_model.output_dir = os.path.abspath(output_dir)
swat_model.out_var_agr = ['sw_mm', 'sw_awc']
swat_model.out_var_hid = ['q_mm/d', 'q_mm/s']
swat_model.results_years = 32
swat_model.spatial_units = 768
swat_model.pdmm = 7
swat_model.verbose = False
swat_model.parallel = True

swat_model.pdmm_param = {'inf_pond':[['POT_FR', 'POT_VOLX'], [0.4, 25] , ['replace','replace'],'hru'],
                      'oil_coff': ['CN2', [45, 66, 77, 80] ,'replace','mgt'], # CN2 values according to soil hydrological group in the hru
                      'corn_cass': ['CN2',[67, 77, 83, 87] ,'replace','mgt'], 
                      'forest': ['CN2', [30, 55, 70, 77] ,'replace','mgt'], 
                      'cha_prot':[['CH_N2','CH_COV2'], [0.15, 0.0,],['replace','replace'],'rte'],
                      'cha_sta':[['CH_S2','CH_COV2'], [-0.1, 0.0],['multiply','replace'],'rte'],
                      'sto_pond': [['WET_FR','WET_NSA','WET_NVOL'], [0.3, 0.5, 1],['replace','replace', 'replace'],'pnd']}


# OPTIMIZATION ROUTINE #################################################################################################

# Optimization settings
seed = 12345       # Seed number (for reproducibility)
n_obj = 2          # Number of objective functions
n_const = 1        # Number of constraints
nparams = 6531     # Number of decision variables
pop_size = 2     # Population size
nmaxgen = 1      # Maximum number of generations (stopping criteria)
opt = 'local'      # Whether the work will be submitted locally or to an HPC
partition = 'serc' # Name of the HPC partition where workers will be deployed

# Step 0: set cluster configuration

if opt == 'local':
    client = Client(processes=False)
elif opt == 'hpc':
    dask.config.set({'distributed.scheduler.allowed-failures': 50})
    cluster = SLURMCluster(cores=1, memory='5G', queue=partition,
                           walltime='00:30:00', processes=1,
                           worker_extra_args=["--lifetime", "25m",
                                              "--lifetime-stagger", "4m"])
    cluster.adapt(minimum_jobs=4, maximum_jobs=pop_size)
    client = Client(cluster)

# Step 1: create optimization problem object
swat_model.n_obj = n_obj
problem = SWATProblem(swat_model, client)

# Step 2: create reference directions
ref_dirs = get_reference_directions("energy", n_obj, pop_size, seed=seed)

# Step 3: create algorithm object
pop_file = open(X_dir)
X_df = pd.read_csv(pop_file, header=None, index_col=False)

X = X_df.to_numpy()
crossover = get_crossover("int_sbx", prob=0.9, eta=10)
mutation = get_mutation("int_pm", eta=20, prob=1/nparams)
algorithm = UNSGA3(ref_dirs=ref_dirs,
                   pop_size=None,
                   sampling=X,
                   crossover=crossover,
                   mutation=mutation,
                   eliminate_duplicates=True,
                   callback=my_callback)

# Step 4: create optimization object
start_time = time.time()
swat_model.clean() # clean any previous unremoved temp. directories
res = minimize(problem,
               algorithm,
               termination=('n_gen', nmaxgen),
               seed=seed,
               verbose=True)
client.close()

# Step 5: report results
print("--- {:.2f} minutes ---".format((time.time() - start_time) / 60))
print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))

# clean temporal folders
swat_model.clean()
