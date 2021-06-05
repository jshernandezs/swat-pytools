#!/usr/bin/env python3
import os
import platform
import time
from pymoo.algorithms.unsga3 import UNSGA3
from pymoo.factory import get_sampling, get_crossover, get_mutation, get_reference_directions
from pymoo.optimize import minimize
from swat_utilities.calibration_factory import SWATConfig, SWATProblem, my_callback
from swat_utilities.custom_problems import set_problem

# create calibration settings object     
swat_cal_config = SWATConfig()
# calibration settings
swat_version = ''
if platform.system() == 'Linux':
    swat_version = 'SWAT_Rev622'
elif platform.system() == 'Darwin':
    swat_version = 'SWAT_Rev622_macOS'

swat_cal_config.model_file = os.path.abspath('../resources/Models/Honeyoy_Model.zip')
swat_cal_config.swat_exec_name = swat_version
swat_cal_config.obs_file = os.path.abspath('../resources/Observed/Honeyoy_cal.csv')
swat_cal_config.out_file = 'watout.dat'
swat_cal_config.out_var = 'FLOWm^3/s'
swat_cal_config.output_dir = os.path.abspath('../output/Problem_1')
swat_cal_config.temp_dir = '/tmp/swat_runs/Problem_1'
swat_cal_config.temp_run_dir = '/tmp/output_swat/Problem_1'
swat_cal_config.verbose = False
swat_cal_config.parallel = True
swat_cal_config.cal_param = {'BIOMIX': [[0, 1], 'replace', 'mgt'],
                             'CN2': [[-0.25, 0.25], 'multiply', 'mgt'],
                             'CANMX': [[0, 100], 'replace', 'hru'],
                             'ESCO': [[0, 1], 'replace', 'hru'],
                             'EPCO': [[0, 1], 'replace', 'hru'],
                             'GW_DELAY': [[0, 500], 'replace', 'gw'],
                             'ALPHA_BF': [[0, 1], 'replace', 'gw'],
                             'GWQMN': [[0, 5000], 'replace', 'gw'],
                             'GW_REVAP': [[0.02, 0.2], 'replace', 'gw'],
                             'REVAPMN': [[0, 1000], 'replace', 'gw'],
                             'RCHRG_DP': [[0, 1], 'replace', 'gw'],
                             'CH_N2': [[0, 0.3], 'replace', 'rte'],
                             'CH_K2': [[0, 500], 'replace', 'rte'],
                             'SOL_AWC': [[-0.25, 0.25], 'multiply', 'sol'],
                             'SURLAG': [[1, 24], 'replace', 'bsn']}

# Optimization routine #######################################################################

# Step 1: create optimization problem object
seed = 12345
swat_cal_config, n_obj = set_problem(swat_cal_config, problem=1)
problem = SWATProblem(swat_cal_config, parallelization=("threads", 7))

# Step 2: create reference directions
pop_size = 100
ref_dirs = get_reference_directions("energy", n_obj, pop_size, seed=seed)

# Step 3: create algorithm object
sampling = get_sampling("real_random")
crossover = get_crossover("real_sbx", prob=0.9, eta=10)
mutation = get_mutation("real_pm", eta=20, prob=1 / 15)
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
               termination=('n_gen', 1000),
               seed=seed,
               verbose=True)

# Step 5: report results
print("--- {:.2f} minutes ---".format((time.time() - start_time) / 60))
print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))

# clean temporal folders
swat_cal_config.clean()
