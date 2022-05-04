#!/usr/bin/env python3
""" Factory of optimization problems using SWAT
"""
import os
import math
import csv
import autograd.numpy as anp
from pymoo.model.problem import Problem
from pymoo.model.individual import Individual
from pymoo.model.population import Population
from pymoo.model.algorithm import filter_optimum
from swat_utilities import swat_config
from multiprocessing import Lock
import pandas as pd
import numpy as np
import subprocess as sp
import datetime as dt


class Hooker:
    """ Stores simulations for each generation/iteration
    """
    simulations = []
    params = []
    obj_f = []
    output_dir = ''

    def add_sim(self, f, param):

        new_param = {key: param[key][0] for key in param.keys()}

        self.obj_f.append(f)
        self.params.append(new_param)

    def print(self, algorithm):
        n_gen = algorithm.n_gen
        if not self.output_dir == '':
            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)

            df2 = pd.DataFrame(self.params)
            df3 = pd.DataFrame(self.obj_f)

            filename2 = 'parameters_gen{:04d}.csv'.format(n_gen)
            filename3 = 'objs_gen{:04d}.csv'.format(n_gen)

            df2.to_csv(self.output_dir + '/' + filename2, header=True, index=False)
            df3.to_csv(self.output_dir + '/' + filename3, header=True, index=False)

            # store optimal solution per generation
            if algorithm.opt is None:
                opt = filter_optimum(algorithm.pop)
            else:
                opt = algorithm.opt

            n_obj = len(algorithm.pop[0].F)

            if isinstance(opt, Population):
                x, f = opt.get("X", "F")
            elif isinstance(opt, Individual):
                x, f = opt.X.reshape(1, -1), opt.F.reshape(1, -1)
            else:
                x, f = None, None

            df4 = pd.DataFrame(x, columns=df2.columns)
            df5 = pd.DataFrame(f, columns=df3.columns[0:n_obj])

            filename4 = 'opt_parameters_gen{:04d}.csv'.format(n_gen)
            filename5 = 'opt_objs_gen{:04d}.csv'.format(n_gen)

            df4.to_csv(self.output_dir + '/' + filename4, header=True, index=False)
            df5.to_csv(self.output_dir + '/' + filename5, header=True, index=False)

        else:
            print('You must define output directory first!')

    @staticmethod
    def clear_hooker():
        Hooker.obj_f = []
        Hooker.params = []

    @staticmethod
    def set_output_dir(out_dir):
        Hooker.output_dir = out_dir


def clean_all(config):
    try:
        temp_dir = config.temp_dir
        temp_run_dir = config.temp_run_dir
        sp.run(['rm', '-rf', temp_dir])
        sp.run(['rm', '-rf', temp_run_dir])
    except sp.CalledProcessError:
        print('temporal folders do not exist')


class SWATConfig:

    def __init__(self):
        self.swat_dir = os.path.abspath('../resources')
        self.model_file = ''
        self.swat_exec_name = 'SWAT_Rev670'
        self.obs_file = ''
        self.out_var_rch = ['FLOW_OUTcms']
        self.out_var_sub = ['SWmm']
        self.cal_param = {}
        self.n_obj = 2
        self.output_dir = os.path.abspath('../output')
        self.temp_dir = '/tmp/swat_runs'
        self.verbose = True
        self.lock = 0
        self.runid = swat_config.RunID()
        self.temp_run_dir = '/tmp/output_swat'

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if not os.path.exists(self.temp_dir):
            os.makedirs(self.temp_dir)

    def clean(self):
        temp_dir = self.temp_dir
        temp_run_dir = self.temp_run_dir
        try:
            sp.run(['rm', '-rf', temp_dir])
            sp.run(['rm', '-rf', temp_run_dir])
        except sp.CalledProcessError:
            print('temporal folders do not exist')


class SWATProblem(Problem):

    def __init__(self, cal_config, **kwargs):
        param = cal_config.cal_param
        n_obj = cal_config.n_obj

        xl = [param[x][0][0] for x in param.keys()]
        xu = [param[x][0][1] for x in param.keys()]

        super().__init__(n_var=len(xl),
                         n_obj=n_obj,
                         xl=anp.array(xl),
                         xu=anp.array(xu),
                         evaluation_of=['F'],
                         elementwise_evaluation=True,
                         **kwargs)

        cal_config.lock = Lock()

        self.cal_config = cal_config

    def _evaluate(self, x, out, *args, **kwargs):
        n_obj = self.n_obj
        cal_config = self.cal_config

        f, f_out, param, simulated = fun(x, cal_config, n_obj)
        out["F"] = anp.array(f)

        # update hooker
        lib = Hooker()
        lib.add_sim(f_out, param)
        lib.set_output_dir(cal_config.output_dir)


def my_callback(algorithm):
    # write hooked values out
    lib = Hooker()
    lib.print(algorithm)
    lib.clear_hooker()
    # clear list of ids
    ids = swat_config.RunID()
    ids.clear_id()


def fun(x, cal_config, n_obj):
    lock = cal_config.lock
    runid = cal_config.runid
    swat_model = swat_config.ModelSetup(cal_config.model_file)
    swat_model.swat_dir = cal_config.swat_dir
    swat_model.temp_dir = cal_config.temp_run_dir
    swat_model.output_dir = cal_config.temp_dir
    swat_model.verbose = cal_config.verbose

    # prepare parameters to change
    param_ref = cal_config.cal_param
    param = {}
    for j, key in enumerate(param_ref.keys()):
        param[key] = [x[j], param_ref[key][1], param_ref[key][2]]

    # swat preparation
    swat_model.swat_exec_name = cal_config.swat_exec_name
    swat_model.param = param

    # swat execution
    simulated = run_single_model(swat_model, cal_config.out_var_rch, cal_config.out_var_sub, lock, runid)

    # objective function computation
    
    f = anp.empty([1, n_obj])
    
    q = simulated[0]
    sm = simulated[1]
      
    f[0,0] = -sum([q_i['AREAkm2'] * q_i['Series'].mean() for q_i in q.values()])/sum([q_i['AREAkm2'] for q_i in q.values()])
    f[0,1] = -sum([sm_i['AREAkm2'] * sm_i['Series'].mean() for sm_i in sm.values()])/sum([sm_i['AREAkm2'] for sm_i in sm.values()])

    f_out = {}
    f_out['obj_1'] = f[0,0]
    f_out['obj_2'] = f[0,1]
    
    return f, f_out, param, simulated

def run_single_model(swat_model, out_var_rch, out_var_sub, lock, runid):
    # assigning ID to swat output folder
    if not isinstance(lock, int):
        lock.acquire()
    ind = runid.get_id()
    swat_model.new_model_name = 'RUN{:04d}'.format(ind)
    # prepare SWAT run
    swat_model.prepare_swat()
    if not isinstance(lock, int):
        lock.release()
    # execute model
    swat_model.run_swat()
    # get output time series of variables of interest
    ts2 = get_sub_output(swat_model, out_var_sub)
    ts1 = get_rch_output(swat_model, out_var_rch)
    simulated = [ts1, ts2]
    # remove temporal output folders
    swat_model.remove_swat()   
    
    return simulated

def get_sub_output(model, out_var):
    
    output_dir = model.output_dir + '/' + model.new_model_name + '/output.sub'
    max_col1 = 35
    max_col2 = 35
    ns = 10
    
    with open(output_dir) as f:
        head = [next(f) for x in range(9)]
        raw_header = head[8]
        part1 = raw_header[:max_col1].split()
        part2 = raw_header[max_col1:]
        
        nvar = int(math.ceil((len(part2)-1)/ns))
        part2 = [part2[int((i-1)*ns):int(i*ns)].strip() for i in range(1,nvar+1)]
        
        header = part1 + part2
        
        reader = csv.reader(f)
        buf = []
        for line in reader: 
            part1 = [line[0][6:11], line[0][11:20], line[0][20:25], line[0][25:35]]
            part1[-1] = '0' + part1[-1]
            part1 = [x.strip() for x in part1]
            part2 = line[0][max_col2:]
            part2 = [float(part2[int((i-1)*ns):int(i*ns)].strip()) for i in range(1,nvar+1)]        
            line_list = part1 + part2
            buf.append(line_list)    
    
    raw = pd.DataFrame(data=buf, columns=header)    
    temp = raw.loc[raw['MON'].astype(float) < 13, :].copy()
    temp['AREAkm2'] = temp['AREAkm2'].astype(float)
    temp2 = raw.loc[raw['MON'].astype(float) > 1000, :].copy()
    years = np.repeat(np.unique(temp2['MON']), 12)
    
    subbasins = np.sort(np.unique(temp['SUB']).astype(int))
    
    output = {}
    for sub in subbasins:
        aux = temp.loc[temp['SUB'] == str(sub), :].copy()
        aux['YEAR'] = years
        aux.index = [dt.datetime(int(year),int(month),1) for year, month in zip(aux['YEAR'], aux['MON'])]
        output[str(sub)] = {'AREAkm2': aux['AREAkm2'][0], 'Series': aux.loc[:, out_var]}    
    
    return output

def get_rch_output(model, out_var):
    
    output_dir = model.output_dir + '/' + model.new_model_name + '/output.rch'
    max_col1 = 37
    max_col2 = 38
    ns = 12
    
    with open(output_dir) as f:
        head = [next(f) for x in range(9)]
        raw_header = head[8]
        part1 = raw_header[:max_col1].split()
        part2 = raw_header[max_col1:]
        
        nvar = int(math.ceil((len(part2)-1)/ns))
        part2 = [part2[int((i-1)*ns):int(i*ns)].strip() for i in range(1,nvar+1)]
        
        header = part1 + part2
        
        reader = csv.reader(f)
        buf = []
        for line in reader: 
            part1 = [line[0][6:11], line[0][11:20], line[0][20:26], line[0][26:38]]
            part1 = [x.strip() for x in part1]
            part2 = line[0][max_col2:]
            part2 = [float(part2[int((i-1)*ns):int(i*ns)].strip()) for i in range(1,nvar+1)]        
            line_list = part1 + part2
            buf.append(line_list)    
    
    raw = pd.DataFrame(data=buf, columns=header)    
    temp = raw.loc[raw['MON'].astype(float) < 13, :].copy()
    temp['AREAkm2'] = temp['AREAkm2'].astype(float)
    temp2 = raw.loc[raw['MON'].astype(float) > 1000, :].copy()
    years = np.repeat(np.unique(temp2['MON']), 12)
    
    subbasins = np.sort(np.unique(temp['RCH']).astype(int))
    
    output = {}
    for sub in subbasins:
        aux = temp.loc[temp['RCH'] == str(sub), :].copy()
        aux['YEAR'] = years
        aux.index = [dt.datetime(int(year),int(month),1) for year, month in zip(aux['YEAR'], aux['MON'])]
        output[str(sub)] = {'AREAkm2': aux['AREAkm2'][0], 'Series': aux.loc[:, out_var]}    
    
    return output