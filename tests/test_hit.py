#!/usr/bin/env python3
import sys
import os
import time

main_dir = os.path.abspath('../main/')
sys.path.insert(1, main_dir)

from swat_utilities.performance import read_obs_file
from swat_utilities.hit import HydrologicIndex as Hi

# read observed time-series
observedFile = '../../resources/Observed/Honeyoy.csv'
observed = read_obs_file(observedFile)
# Compute index of interest
start_time = time.time()
hit_obs = Hi(observed, drainage_area=1061.31)
index = 'MAG5'
ind = hit_obs.get_index(index)
print("--- Index {:s}: {:.5f} seconds ---".format(index, (time.time() - start_time)))
# Compute all indices
start_time = time.time()
hit_obs = Hi(observed, drainage_area=1061.31)
hit = hit_obs.get_all_hit_indices()
print("--- All HIT indices: {:.5f} seconds ---".format((time.time() - start_time)))
start_time = time.time()
mag7 = hit_obs.get_mag_seven()
print("--- Magnificent seven indices: {:.5f} seconds ---".format((time.time() - start_time)))
# Compute list of indices
start_time = time.time()
hit_obs = Hi(observed, drainage_area=1061.31)
list_ind = ['dh17', 'mag5', 'fl1', 'ma41', 'ma43', 'ma45']
hit_list = hit_obs.get_list_indices(list_ind)
print("--- List of indices: {:.5f} seconds ---".format((time.time() - start_time)))
