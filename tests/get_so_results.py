#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np

from swat_utilities.post_processing import Optimal
from swat_utilities.performance import Metric, read_obs_file

metrics = ['KGE', 'KGEp', 'NSE', 'IoA', 'R2', 'ksR2', 'R4MS4E']
opt_transform = [True, True, True, True, True, True, False]
transforms = ['none', 'sqrt', 'log', 'rel', 'inverse']
observed_file = os.path.abspath('../resources/Observed/Honeyoy.csv')  # full date

c = 0
param_all = []
sim_all = []
of_all = []
names = []
pm_all = []
pm_names = []
metric_names = ['NSE', 'RMSE', 'KGE', 'KGEp', 'RSR', 'IoA', 'R2', 'koR2', 'ksR2', 'R4MS4E', 'PBIAS']

for i, metric in enumerate(metrics):
    for transform in transforms:
        if opt_transform[i] or (not opt_transform[i] and transform == 'none'):
            cal_dir = '../output/Single_Obj/test_{:s}_{:s}'.format(metric, transform)
            names.append('{:s}_{:s}'.format(metric, transform))
            # create object with optimal solutions
            swat_optimal = Optimal(cal_dir)
            # obtain optimal solutions
            best_of, best_param, best_sim = swat_optimal.get_best()
            # set dates range for simulations based on observed range
            obs = read_obs_file(observed_file)
            best_sim = best_sim[obs.index[0]:obs.index[-1]]
            # append results to list of lists or data frames
            param_all.append(best_param)
            of_all.append(best_of.values)
            if c == 0:
                sim_all = best_sim
            else:
                sim_all = pd.merge(sim_all, best_sim, how='inner', left_index=True, right_index=True)
            # compute metrics for each optimal solution
            simulation = Metric(observed_file, best_sim)
            pm = ()
            pm_names = []
            k = 0
            for tr in transforms:
                aux = simulation.get_all(tr)
                aux2 = ['{:s}_{:s}'.format(x, tr) for x in metric_names]
                if k != 0:
                    aux = aux[: len(aux) - 2]
                    aux2 = aux2[: len(aux2) - 2]
                pm = pm + aux
                pm_names = pm_names + aux2
                k += 1
            pm_all.append(list(pm))
            c += 1
        else:
            continue

# convert aggregated results to data frames
metrics_all = pd.DataFrame(pm_all, index=names, columns=pm_names)
param_all = pd.concat(param_all, axis=0)
param_all.index = names
of_all = pd.DataFrame(np.array(of_all).reshape(-1, 1))
of_all.index = names
sim_all.columns = names

# write results to .csv files
output_dir = os.path.abspath('../optimal')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

metrics_all.to_csv(output_dir + '/' + 'performance.csv')
param_all.to_csv(output_dir + '/' + 'parameters.csv')
of_all.to_csv(output_dir + '/' + 'objectives.csv', header=False)
sim_all.to_csv(output_dir + '/' + 'simulations.csv')
