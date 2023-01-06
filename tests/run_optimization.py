#!/usr/bin/env python3
"""
Generic optimization routine using monthly streamflow and soil water outputs
from SWAT

"""

import os
import time
from pymoo.algorithms.moo.unsga3 import UNSGA3
from pymoo.factory import get_sampling, get_crossover, get_mutation, get_reference_directions
from pymoo.optimize import minimize
from swat_utilities.optimization_factory import SWATConfig, SWATProblem, my_callback
import dask
from dask_jobqueue import SLURMCluster
from dask.distributed import Client

# MODELING SETTINGS ####################################################################################################

output_dir = '../resources/swat_output/Test_AnaMaria'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 1) create calibration settings object
swat_model = SWATConfig()

# 2) General settings
swat_model.model_file = os.path.abspath('../resources/Models/Test.zip')
swat_model.agr_treshold = os.path.abspath('../resources/csv_files/treshold_agr_ecdf_30.csv')
swat_model.hid_treshold = os.path.abspath('../resources/csv_files/treshold_hid_ecdf_30.csv')
swat_model.swat_exec_name = 'SWAT_Rev681'
swat_model.temp_dir = '/tmp/swat_runs'
swat_model.temp_run_dir = '/tmp/output_swat'
swat_model.output_dir = os.path.abspath(output_dir)
swat_model.out_var_sub = ['SWmm', 'severity']
swat_model.out_var_rch = ['FLOW_OUTcms', 'severity']
swat_model.verbose = False
swat_model.parallel = True
swat_model.cal_param = {'POT_FR': [[0, 0.3], 'replace', 'hru'],
                        'POT_VOLX': [[0, 300], 'replace', 'hru'],
                        'CN2': [[-0.4, 0.4], 'multiply', 'mgt'],
                        'CHW2': [[-0.2, 0.2], 'multiply', 'rte'],
                        'CH_S2': [[-0.2, 0.2], 'multiply', 'rte'],
                        'CH_N2': [[-0.2, 0.2], 'multiply', 'rte'],
                        'CH_S1': [[-0.2, 0.2], 'multiply', 'sub'],
                        'CH_W1': [[-0.2, 0.2], 'multiply', 'sub'],
                        'CH_N1': [[-0.2, 0.2], 'multiply', 'sub'],
                        'PND_FR': [[0, 0.3], 'replace', 'pnd'],
                        'PND_PSA': [[0, 5], 'replace', 'pnd'],
                        'PND_PVOL': [[0, 5], 'replace', 'pnd'],
                        'WET_FR': [[0, 0.3], 'replace', 'pnd'],
                        'WET_NSA': [[0, 5], 'replace', 'pnd'],
                        'WET_NVOL': [[0, 5], 'replace', 'pnd']}

# OPTIMIZATION ROUTINE #################################################################################################

# Optimization settings
seed = 12345       # Seed number (for reproducibility)
n_obj = 2          # Number of objective functions
nparams = 15       # Number of decision variables
pop_size = 7       # Population size
nmaxgen = 2        # Maximum number of generations (stopping criteria)
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
sampling = get_sampling("real_random")
crossover = get_crossover("real_sbx", prob=0.9, eta=10)
mutation = get_mutation("real_pm", eta=20, prob=1/nparams)
algorithm = UNSGA3(ref_dirs=ref_dirs,
                   pop_size=None,
                   sampling=sampling,
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
