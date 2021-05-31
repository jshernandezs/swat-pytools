#!/usr/bin/env python3
"""
Routine for selecting performance metrics based on optimization results and ERHIs replication
"""

import sys
import os
import pandas as pd
import numpy as np
import datetime as dt

mainDIR = os.path.abspath('../main/')
sys.path.insert(1, mainDIR)

from swat_utilities.performance import read_obs_file
from swat_utilities.hit import HydrologicIndex as Hi

  
def main(file_dir, obs_file, dict_pm, n_max=6):
        
    flow = pd.read_csv(file_dir, index_col=0)
    dates = [dt.datetime.strptime(t, '%Y-%m-%d') for t in flow.index]
    flow.index = dates
    
    obs = read_obs_file(obs_file)
    obs_hit = Hi(obs)
    index_obs_dict = obs_hit.get_all_hit_indices()
    index_obs_dict.update({'MAG': obs_hit.get_mag_seven()})
    index_obs = dict2list(index_obs_dict)
    ind_names = get_ind_names(index_obs_dict)
    
    flow = flow[obs.index[0]:obs.index[-1]]
    
    pm_lab = []
    errors = {}
    for item in dict_pm.items():
        pm = item[0]
        for transform in item[1]:
            name = '{:s}_{:s}'.format(pm, transform)
            pm_lab.append('{:s} ({:s})'.format(pm, transform))
            sim = pd.DataFrame(flow[name])
            sim_hit = Hi(sim, obs_file=obs_file)
            error = array_error_all_indices(sim_hit, index_obs)
            errors[name] = error
    
    error_df = pd.DataFrame(errors, index=ind_names)
    list_pm, bad_ind, good_ind = select_pm(error_df, n_max=n_max)
    
    return (list_pm, bad_ind, good_ind, error_df)

def filter_pm(dict_pm, file):
    mat_perf = pd.read_csv(file, index_col=0)
    list_pm = []
    for item in dict_pm.items():
        pm = item[0]
        for transform in item[1]:
            name = '{:s}_{:s}'.format(pm, transform)
            list_pm.append(name)
    
    mat = mat_perf.copy()
    mat = mat.loc[list_pm, list_pm]
    
    new_dict = {}
    surrogates = []
    
    for item in dict_pm.items():
        pm = item[0]
        new_t = []
        iflag = 0
        for transform in item[1]:
            name = '{:s}_{:s}'.format(pm, transform)
            temp = mat[name].copy()
            if pm not in ['RM4MS4E']:
                best_pm = temp[temp == temp.max()].index.to_list()[0]
            else:
                best_pm = temp[temp == temp.min()].index.to_list()[0]
            
            if name == best_pm:
                iflag = 1
                new_t.append(transform)
            else:
                surrogates.append([name, best_pm])
        
        if iflag == 1:
            new_dict[pm] = new_t
            
    replaced = [x[0] for x in surrogates]
    replacer = [x[1] for x in surrogates]
    
    restore = []
    for item in replaced:
        if item in replacer:
            restore.append(item)
    
    
    for item in restore:
        pm = item.split('_')[0]
        transform = item.split('_')[1]
        if pm in new_dict:
            l0 = new_dict[pm]
            l0.append(transform)
            new_dict[pm] = l0
        else:
            new_dict[pm] = [transform]
    
    return new_dict, surrogates

def select_pm(df, thresh=30.0, n_max=6):
    """
    
    Parameters
    ----------
    df : Pandas DataFrame
        Dataframe containing relative errors for each performance metrics (columns) and hydrologic index (rows)
    thresh : float, optional
        Relative error upper threshold - absolute value (%) to consider an error as acceptable. The default is 30.
    n_max : int, optional
        Maximum number of performance metrics to be selected. The default is 6.

    Returns
    -------
    list_pm : list
        List containing name of selected performance metrics
    bad_ind: list
        List of not well represented ERHI

    """ 
    i = 0
    mat = abs(df.copy())
    mat.iloc[mat <= thresh] = 1
    mat.iloc[mat > thresh] = 0
    list_pm = [[]]
    good_ind = [[]]
    mat = [mat]
    iflag = 0
    while (i < n_max) and (iflag == 0):
        new_m = []
        new_list = []
        new_good = []
        for j, m in enumerate(mat):
            list0 = list_pm[j]
            good0 = good_ind[j]
            n_well = m.sum(axis=0)
            best_pm = n_well[n_well == n_well.max()].index.to_list()
            
            if n_well.max() == 0:
                iflag = 1
                break
            
            for pm in best_pm:
                list1 = list0 + [[pm, n_well.max()]]
                error_ref = m[pm]
                
                good1 = good0 + [error_ref[error_ref == 1].index.to_list()]
                
                m1 = m.drop(pm, axis=1)
                m1 = m1.drop(error_ref[error_ref == 1].index, axis=0)   
                
                new_m.append(m1)
                new_list.append(list1)
                new_good.append(good1)
        
        if iflag == 0:        
            list_pm = new_list
            good_ind = new_good
            mat = new_m
            i += 1
    
    bad_ind = [x.index.to_list() for x in mat]
    
    return (list_pm, bad_ind, good_ind)

def dict2list(d):
    i = 0
    for k, inds in d.items():
        if i == 0:
            l = inds
        else:
            l = l + inds
        i += 1
    return(l)

def array_error_all_indices(sim_hit, index_obs):
    
    eps = np.finfo(float).eps   
    
    index_sim_dict = sim_hit.get_all_hit_indices()
    index_sim_dict.update({'MAG': sim_hit.get_mag_seven()})    
    index_sim = dict2list(index_sim_dict)
    
    error = [((s+eps) - (o+eps)) / (o+eps) * 100 for o, s in zip(index_obs, index_sim)]    
    return(error)

def get_ind_names(d):
    d_names = []
    for k, ind in d.items():
        names = ['{:s}{:d}'.format(k.lower(), i) for i in range(1, len(ind)+1)]
        d_names = d_names + names
    return(d_names)

if __name__ == "__main__":
    
    file_dir = '../../optimal/simulations.csv'
    performance_file = '../../optimal/performance.csv'
    obs_file = '../../resources/Observed/Honeyoy.csv'
    dict_pm = {'KGE': ['none', 'sqrt', 'log', 'inverse'],
               'KGEp': ['none', 'sqrt', 'log', 'inverse'],
               'NSE': ['none', 'sqrt', 'log', 'rel', 'inverse'],
               'R2': ['none', 'sqrt', 'log', 'inverse'],
               'IoA': ['none', 'sqrt', 'log', 'rel', 'inverse'],
               'R4MS4E': ['none']}
    
    dict_pm_f, surrogates  = filter_pm(dict_pm, performance_file)
    
    list_pm, bad_ind, good_ind, error_df = main(file_dir, obs_file, dict_pm, n_max=6)
    list_pm_f, bad_ind_f, good_ind_f, error_df_f = main(file_dir, obs_file, dict_pm_f, n_max=6)
    