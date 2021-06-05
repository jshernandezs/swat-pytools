#!/usr/bin/env python3
"""
Run multiple SWAT models from a list of given parameters

"""
import os
import platform
import copy as cp
import pickle
from multiprocessing.pool import ThreadPool
from swat_utilities.swat_config import ModelSetup
from swat_utilities.post_processing import Optimal
   

def main(cal_dirs, model_file_path, params_template, processes=4):
    
    if platform.system() == 'Linux':
        swat_version = 'SWAT_Rev622'
    elif platform.system() == 'Darwin':
        swat_version = 'SWAT_Rev622_macOS'
    else:
        swat_version = None
        print('Windows version is not supported yet')
    
    swat_model = ModelSetup(model_file_path)
    swat_model.swat_exec_name = swat_version
    swat_model.verbose_swat = False
    
    ret = []
    for cal_dir in cal_dirs:
        temp = run_batch(swat_model, cal_dir, params_template, processes)
        ret.append(temp)
    
    return ret

def get_sim(model, params, name):
    
    m = cp.deepcopy(model)
    m.new_model_name = name
    m.param = [params]
    m.prepare_swat()
    m.run_swat()
    m.read_output()
    m.simulated('watout.dat', plot=False)
    sim = m.sim_ts['FLOWm^3/s']
    
    m.remove_swat()
    
    return sim

def run_batch(swat_model, cal_dir, params_template, processes):
    
    swat_optimal = Optimal(cal_dir)
    params_set = swat_optimal.get_params()
    
    list_new_params = []
    for i, pars in params_set.iterrows():
        new_params = params_template.copy()
        for p_name, item in new_params.items():
            new_params[p_name][0] = float(pars[p_name])
        list_new_params.append(cp.deepcopy(new_params))
    
    inputs = [(swat_model, new_param, 'RUN_{:04d}'.format(i)) for i, new_param in enumerate(list_new_params)]    
    
    with ThreadPool(processes) as pool:
        sim_list = pool.starmap(get_sim, inputs)    
    
    i = 0
    for sim in sim_list:
        sim.columns = [str(i)]
        if i == 0:
            df = sim
        else:                
            df = df.join(sim)
        i += 1
        
    return df 
    

if __name__ == '__main__':
    
    cal_dirs = ['../output/Performance_based', '../output/Signature_based']
    model_file_path = os.path.abspath('../resources/Models/Honeyoy_Model_val.zip')
    
    params_template = {'BIOMIX': [0.22, 'replace', 'mgt'],
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
    
    simulations = main(cal_dirs, model_file_path, params_template, 8)
    
    with open('../output/validation.pickle', 'wb') as f:
        pickle.dump(simulations, f, pickle.HIGHEST_PROTOCOL)    
