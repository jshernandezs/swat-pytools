#!/usr/bin/env python3
"""
Wrapper for executing SWAT from Matlab
"""

import sys
from swat_utilities.swat_config import ModelSetup
from swat_utilities.performance import Metric

def main(x, model_name):

    model_file_path = '../../../resources/Models/Honeyoy_Model_val.zip';
    
    swat_model = ModelSetup(model_file_path)
    
    swat_model.swat_dir = '../../../resources'
    swat_model.swat_exec_name = 'SWAT_Rev622'
    swat_model.new_model_name = model_name
    swat_model.verbose_swat = False
    swat_model.verbose = False
    
    observed_file = '../../../resources/Observed/Honeyoy_val.csv'
    outputfile = 'watout.dat'
    
    swat_model.param = {
                        'BIOMIX': [x[0], 'replace', 'mgt'],
                        'CN2': [x[1], 'multiply', 'mgt'],
                        'CANMX': [x[2], 'replace', 'hru'],
                        'ESCO': [x[3], 'replace', 'hru'],
                        'EPCO': [x[4], 'replace', 'hru'],
                        'GW_DELAY': [x[5], 'replace', 'gw'],
                        'ALPHA_BF': [x[6], 'replace', 'gw'],
                        'GWQMN': [x[7], 'replace', 'gw'],
                        'GW_REVAP': [x[8], 'replace', 'gw'],
                        'REVAPMN': [x[9], 'replace', 'gw'],
                        'RCHRG_DP': [x[10], 'replace', 'gw'],
                        'CH_N2': [x[11], 'replace', 'rte'],
                        'CH_K2': [x[12], 'replace', 'rte'],
                        'SOL_AWC': [x[13], 'multiply', 'sol'],
                        'SURLAG': [x[14], 'replace', 'bsn']
                        }
    
    swat_model.prepare_swat()
    swat_model.run_swat()
    
    swat_model.read_output()
    swat_model.simulated(outputfile, plot=False);
    simulation = Metric(observed_file, swat_model.sim_ts['FLOWm^3/s'])
    
    simulation.series.to_csv(model_name + '.dat', header=False, index=False)
    swat_model.remove_swat()    

if __name__ == '__main__':
    
    x_raw = sys.argv[1]
    name = sys.argv[2]
    
    n = len(sys.argv[1])
    x = sys.argv[1][1:n-1].split(',')
    x = [float(i) for i in x]
    
    main(x, name)