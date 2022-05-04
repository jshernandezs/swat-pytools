#!/usr/bin/env python3
"""
Generic optimization routine using monthly streamflow and soil water outputs
from SWAT

@author: J. Sebastian Hernandez-Suarez
"""

import os
import time
import sys
from pymoo.algorithms.unsga3 import UNSGA3
from pymoo.factory import get_sampling, get_crossover, get_mutation, get_reference_directions
from pymoo.optimize import minimize

main_dir = '/home/jshs/Repositories/swat-pytools/src/' # <-- swat_utilities directory 
sys.path.insert(1, main_dir)

from swat_utilities.optimization_factory import SWATConfig, SWATProblem, my_callback

# MODELING SETTINGS ####################################################################################################

output_dir = '../resources/swat_output/Test_AnaMaria'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 1) create calibration settings object     
swat_model = SWATConfig()

# 2) General settings
swat_model.model_file = os.path.abspath('../resources/Models/Test.zip')
swat_model.swat_exec_name = 'SWAT_Rev681'
swat_model.output_dir = os.path.abspath(output_dir)
swat_model.out_var_rch = ['FLOW_OUTcms']
swat_model.out_var_sub = ['SWmm']
swat_model.verbose = False
swat_model.parallel = True
swat_model.cal_param = {'CANMX': [[0, 100], 'replace', 'hru'],
                        'ESCO': [[0, 1], 'replace', 'hru'],
                        'EPCO': [[0, 1], 'replace', 'hru']}

# OPTIMIZATION ROUTINE #################################################################################################

# Optimization settings
seed = 12345    # Seed number (for reproducibility)
n_obj = 2       # Number of objective functions
nparams = 3     # Number of decision variables
nprocesses = 4  # Number of processes to run in parallel
pop_size = 4    # Population size
nmaxgen = 2     # Maximum number of generations (stopping criteria)

# Step 1: create optimization problem object
swat_model.n_obj = n_obj
problem = SWATProblem(swat_model, parallelization=("threads", nprocesses))

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
res = minimize(problem,
               algorithm,
               termination=('n_gen', nmaxgen),
               seed=seed,
               verbose=True)

# Step 5: report results
print("--- {:.2f} minutes ---".format((time.time() - start_time) / 60))
print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))

# clean temporal folders
swat_model.clean()