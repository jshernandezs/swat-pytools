#!/usr/bin/env python3
"""
Summary of single-objective calibrations in terms of ERHIs
"""

import os
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# import seaborn as sns

from swat_utilities.performance import read_obs_file
from swat_utilities.hit import HydrologicIndex as Hi
from swat_utilities.post_processing import heatmap


file_dir = '../optimal/simulations.csv'
obs_file = '../resources/Observed/Honeyoy.csv'
dict_pm = {'KGE': ['none', 'sqrt', 'log', 'inverse'],
           'KGEp': ['none', 'sqrt', 'log', 'inverse'],
           'NSE': ['none', 'sqrt', 'log', 'rel', 'inverse'],
           'R2': ['none', 'sqrt', 'log', 'inverse'],
           'IoA': ['none', 'sqrt', 'log', 'rel', 'inverse'],
           'R4MS4E': ['none']}

# list_ind = [['ma12', 'ma13', 'ma14', 'ma15', 'ma16', 'ma17', 'ma18', 'ma19', 'ma20', 'ma21', 'ma22', 'ma23'],
#             ['dl1', 'dl2', 'dl3', 'dl4', 'dl5', 'dh1', 'dh2', 'dh3', 'dh4', 'dh5', 'ml17'], ['tl1', 'th1',
#             'fl1', 'dl16', 'fh1', 'dh15', 'ra1', 'ra3', 'ra8'], 
#             ['mag1', 'mag2', 'mag3', 'mag4', 'mag5', 'mag6', 'mag7']]

# groups = ['IHA1', 'IHA2', 'IHA3_4_5', 'MAG7']
# titles = ['Magnitude (monthly)',
#           'Magnitude and duration (annual extremes)',
#           'Duration, frequency, rate of change, and timing',
#           'Magnificent seven']

list_ind = [['ma{:d}'.format(x) for x in range(1,46)],
            ['ml{:d}'.format(x) for x in range(1,23)],
            ['mh{:d}'.format(x) for x in range(1,28)],
            ['dl{:d}'.format(x) for x in range(1,21)],
            ['dh{:d}'.format(x) for x in range(1,25)],
            ['fl{:d}'.format(x) for x in range(1,4)],
            ['fh{:d}'.format(x) for x in range(1,12)],
            ['ta{:d}'.format(x) for x in range(1,4)],
            ['tl{:d}'.format(x) for x in range(1,5)],
            ['th{:d}'.format(x) for x in range(1,4)],
            ['ra{:d}'.format(x) for x in range(1,10)],
            ['mag{:d}'.format(x) for x in range(1,8)]]

groups = ['MA', 'ML', 'MH', 'DL', 'DH', 'FL', 'FH', 'TA', 'TL', 'TH', 'RA', 'MAG']
titles = ['a) MA', 'b) ML', 'c) MH', 'd) DL', 'e) DH', 'f) FL', 'g) FH', 'h) TA', 'i) TL', 'j) TH', 'k) RA', 'l) MAG']

pm_lab_f = ['NSE', 'KGE', 'KGE\'', 'IoA', 'R4MS4E', 
            '$NSE_{sqrt}$', '$KGE_{sqrt}$', '$KGE\'_{sqrt}$',
            '$IoA_{sqrt}$', '$KGE_{log}$', '$KGE\'_{log}$', '$IoA_{log}$', 
            '$NSE_{log}$', '$NSE_{inv}$', '$KGE_{inv}$', '$KGE\'_{inv}$', '$IoA_{inv}$',
            '$R^2$', '$R^2_{sqrt}$', '$R^2_{log}$', '$R^2_{inv}$',
            '$NSE_{rel}$', '$IoA_{rel}$']

name_f = ['NSE_none', 'KGE_none', 'KGEp_none', 'IoA_none', 'R4MS4E_none', 
          'NSE_sqrt', 'KGE_sqrt', 'KGEp_sqrt', 
          'IoA_sqrt', 'KGE_log', 'KGEp_log', 'IoA_log',
          'NSE_log', 'NSE_inverse', 'KGE_inverse', 'KGEp_inverse', 'IoA_inverse',
          'R2_none', 'R2_sqrt', 'R2_log', 'R2_inverse',
          'NSE_rel', 'IoA_rel']

output_dir = '../figures'

if not os.path.exists(output_dir):
        os.makedirs(output_dir)

flow = pd.read_csv(file_dir, index_col=0)
dates = [dt.datetime.strptime(t, '%Y-%m-%d') for t in flow.index]
flow.index = dates

obs = read_obs_file(obs_file)
obs_hit = Hi(obs)

flow = flow[obs.index[0]:obs.index[-1]]

all_errors = {}
eps = np.finfo(float).eps

fig = plt.figure(figsize=(12,25))
G = gridspec.GridSpec(4, 8, figure=fig, width_ratios=[3,11,3,4,3,9,7,1], hspace=0.2, wspace=0.1)
axs = [plt.subplot(G[0,:7]),
       plt.subplot(G[1,:3]), plt.subplot(G[1,3:7]),
       plt.subplot(G[2,:3]), plt.subplot(G[2,3:7]), 
       plt.subplot(G[3,0]), plt.subplot(G[3,1]), plt.subplot(G[3,2]),
       plt.subplot(G[3,3]), plt.subplot(G[3,4]), plt.subplot(G[3,5]), plt.subplot(G[3,6]),
       plt.subplot(G[1:3,7])]

c = 0

for t, group, inds in zip(titles, groups, list_ind):
    pm_lab = []
    errors = {}
    for item in dict_pm.items():
        pm = item[0]
        for transform in item[1]:
            name = '{:s}_{:s}'.format(pm, transform)
            pm_lab.append('{:s} ({:s})'.format(pm, transform))
            sim = pd.DataFrame(flow[name])
            sim_hit = Hi(sim, obs_file=obs_file)
            index_obs_dict = obs_hit.get_list_indices(inds)
            index_sim_dict = sim_hit.get_list_indices(inds)
            error = [((s+eps) - (o+eps)) / (o+eps) * 100 
                     for o, s in zip(list(index_obs_dict.values()), list(index_sim_dict.values()))]
            errors[name] = error
    
    errors = pd.DataFrame(errors, index=list(index_obs_dict.keys()))
    inds_plot = [str(i) for i in range(1, len(inds) + 1)]
    errors = errors[name_f]
    
    ax = axs[c]
    
    if c in [0,1,3,5]:
        im = heatmap(errors.to_numpy().transpose(), pm_lab_f, inds_plot, ax=ax,
                           cmap="RdYlBu_r", title = t,
                           vmin=-40, vmax=40, aspect = 'auto')
    else:
        im = heatmap(errors.to_numpy().transpose(), [], inds_plot, ax=ax,
                           cmap="RdYlBu_r", title = t,
                           vmin=-40, vmax=40, aspect = 'auto')
    
    all_errors[group] = errors
    c += 1
    
cb = fig.colorbar(im, cax=axs[c], extend='both')
cb.ax.set_ylabel("Relative error [%]", rotation=-90, va="bottom", fontsize=18)

plt.savefig(output_dir + '/Heatmaps.tif', dpi=600, bbox_inches='tight', pil_kwargs={"compression": "tiff_lzw"})

# Summary per category
ind_numbers = {}
for cat, errors in all_errors.items():
    temp = errors.copy()
    temp[np.abs(errors)<=30] = 1
    temp[np.abs(errors)>30] = 0
    ind_numbers[cat] = temp.sum(axis=0)

ind_perf = pd.DataFrame(ind_numbers)
    
# # Clustering analysis for heatmaps
# for item in all_errors.values():
#     errors = item.copy()
#     errors[errors<-30] = -30
#     errors[errors>30] = 30
#     g = sns.clustermap(errors, cmap='RdYlBu_r', method='ward', yticklabels=1)
