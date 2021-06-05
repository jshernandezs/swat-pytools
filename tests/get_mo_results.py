#!/usr/bin/env python3
import os

from swat_utilities.post_processing import Optimal

obs_file = os.path.abspath('../resources/Observed/Honeyoy.csv')
master_dir = '../output/'

dir_list = os.listdir(master_dir)

for folder in dir_list:
    cal_dir = master_dir + folder
    output_dir = os.path.abspath('../optimal/' + folder)
    swat_optimal = Optimal(cal_dir)

    best_of, best_param, best_sim = swat_optimal.get_best()
    swat_optimal.write_best(output_dir, obs_file)
