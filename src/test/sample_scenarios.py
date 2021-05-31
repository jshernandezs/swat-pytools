#!/usr/bin/env python3
"""
Routine for creating folders with sample files for all possible scenarios
"""

import sys
import os
import time
import multiprocessing
import platform
import subprocess as sp
import pandas as pd

main_dir = os.path.abspath('../main/')
sys.path.insert(1, main_dir)

import swat_utilities.rotation_factory as builder
from swat_utilities.swat_config import ModelSetup


def single_run(output_dir, comb, start_fs, update_params, subbasins, swat_version, model_file_path):
    folder_name = '{:d}{:02d}{:d}{:d}{:d}{:d}'.format(*comb)
    rotation, params, filter_strip = builder.set_rot(*comb, start_fs)
    
    update_params.update(params)
    
    # build SWAT model object
    swat_model = ModelSetup(model_file_path)
    swat_model.swat_exec_name = swat_version
    swat_model.param = [update_params]
    swat_model.drainmod = True
    swat_model.mgt_table = [rotation]
    swat_model.ops_table = [filter_strip]
    swat_model.subbasins = [subbasins]
    swat_model.subbasins_mgt = [subbasins]
    swat_model.subbasins_ops = [subbasins]
    swat_model.new_model_name = folder_name
    
    # run the SWAT model
    swat_model.prepare_swat()
    swat_model.run_swat()  

    # write outputs
    new_dir = output_dir + '/' + folder_name
    output_file = swat_model.output_dir +  '/' + folder_name + '/output.hru'

    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    with open(os.path.abspath(new_dir + '/rotation.dat'), "w") as f:
        f.writelines([x + '\n' for x in rotation.op_sched])

    with open(os.path.abspath(new_dir + '/bmp.dat'), "w") as f:
        f.writelines([x + '\n' for x in filter_strip.schedule])

    with open(os.path.abspath(new_dir + '/parameters.dat'), "w") as f:
        output = ['{:s}: {:s}\n'.format(k, str(v)) for k, v in update_params.items()]
        f.writelines(output)
    
    sp.run(['mv', output_file, new_dir], check=True)
    
    # remove SWAT model folder
    swat_model.remove_swat()


def run():
    
    processes = 20
    output_dir = '/home/sebwaymaster/Desktop/Saginaw_Project/saginaw_scenarios'
    subbasins_filename = '../../resources/csv_files/cropland.csv'
    
    model_file_path = os.path.abspath('../../resources/Models/Saginaw_Model_WQ.zip')
    if platform.system() == 'Linux':
        swat_version = 'SWAT_Rev670'
    elif platform.system() == 'Darwin':
        swat_version = 'SWAT_Rev622_macOS'
    else:
        swat_version = None
        print('Windows version is not supported yet')
        
    new_params = {'NBYR': [12, 'replace', 'cio'],
                  'IYR': [1998, 'replace', 'cio'],
                  'IDAL': [365, 'replace', 'cio'],
                  'NYSKIP': [2, 'replace', 'cio'],
                  'IPRINT': [2, 'replace', 'cio'],
                  'ICALEN': [0, 'replace', 'cio']}

    opt1 = [x for x in range(1, 3)]  # Rotation
    opt2 = [x for x in range(1, 11)]  # System
    opt3 = [x for x in range(1, 6)]  # Application rate
    opt4 = [x for x in range(1, 6)]  # Cover crop
    opt5 = [x for x in range(1, 3)]  # Tile drainage option
    opt6 = [x for x in range(1, 3)]  # Filter strip option

    start_fs = {'MONTH': 1, 'DAY': 1, 'IYEAR': 1998}
    combs = [[a, b, c, d, e, f] for a in opt1 for b in opt2 for c in opt3 for d in opt4 for e in opt5 for f in opt6]
    subbasins = pd.read_csv(subbasins_filename, header=None).iloc[:, 0].to_list()

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with multiprocessing.Pool(processes) as pool:
        results = [pool.apply_async(single_run, args=(output_dir, x, start_fs, new_params, subbasins, swat_version, 
                                                      model_file_path)) for x in combs]
        for r in results:
            r.get()


if __name__ == '__main__':
    multiprocessing.freeze_support()
    start_time = time.time()
    run()
    print("--- Scenarios simulation time was {:.2f} minutes ---".format((time.time() - start_time) / 60))
