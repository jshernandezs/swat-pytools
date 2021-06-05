#!/usr/bin/env python3
import os
import platform
from pymoo.algorithms.so_genetic_algorithm import GA
from pymoo.factory import get_crossover, get_mutation
from pymoo.optimize import minimize
from swat_utilities.calibration_factory import SWATConfig, SWATProblem, my_callback

# create calibration settings object     
swat_cal_config = SWATConfig()
# calibration settings
if platform.system() == 'Linux':
    swat_version = 'SWAT_Rev622'
elif platform.system() == 'Darwin':
    swat_version = 'SWAT_Rev622_macOS'
else:
    swat_version = None
    print('Windows version is not supported yet')

swat_cal_config.model_file = os.path.abspath('../resources/Models/Honeyoy_Model.zip')
swat_cal_config.swat_exec_name = swat_version
swat_cal_config.obs_file = os.path.abspath('../resources/Observed/Honeyoy_cal.csv')
swat_cal_config.out_file = 'watout.dat'
swat_cal_config.out_var = 'FLOWm^3/s'
swat_cal_config.output_dir = os.path.abspath('../output/test_KGE')
swat_cal_config.temp_dir = '/tmp/swat_runs/KGE_none'
swat_cal_config.temp_run_dir = '/tmp/output_swat/KGE_none'
swat_cal_config.obj_f = {'KGE': 'none'}
swat_cal_config.verbose = False
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
problem = SWATProblem(swat_cal_config, parallelization=("threads", 7))
# Step 2: create algorithm object
crossover = get_crossover("real_sbx", prob=0.9, eta=10)
mutation = get_mutation("real_pm", eta=20, prob=1 / 15)
algorithm = GA(pop_size=100,
               crossover=crossover,
               mutation=mutation,
               eliminate_duplicates=True,
               callback=my_callback)
# Step 3: create optimization object
res = minimize(problem,
               algorithm,
               termination=('n_gen', 250),
               seed=1,
               verbose=True)
# Step 4: report results
print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))

# clean temporal folders
swat_cal_config.clean()
