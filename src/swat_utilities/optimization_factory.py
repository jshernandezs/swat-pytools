#!/usr/bin/env python3
""" Factory of optimization problems using SWAT
"""
import os
import math
import csv
import autograd.numpy as anp
from pymoo.core.problem import Problem
from pymoo.core.individual import Individual
from pymoo.core.population import Population
from pymoo.util.optimum import filter_optimum
from swat_utilities import swat_config
import pandas as pd
import numpy as np
import subprocess as sp
import datetime as dt
from datetime import datetime


class Hooker:
    """ Stores simulations for each generation/iteration
    """
    simulations = []
    params = []
    obj_f = []
    output_dir = ''

    def add_sim(self, f, param):

        for f_i, param_i in zip(f, param):
            new_param = {key: param_i[key][0] for key in param_i.keys()}

            self.obj_f.append(f_i)
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
        self.agr_treshold = os.path.abspath('../resources/csv_files/treshold_agr_ecdf_30.csv')
        self.hid_treshold = os.path.abspath('../resources/csv_files/treshold_hid_ecdf_30.csv')
        self.swat_exec_name = 'SWAT_Rev670'
        self.obs_file = ''
        self.out_var_rch = ['FLOW_OUTcms', 'severity']
        self.out_var_sub = ['SWmm', 'severity']
        self.cal_param = {}
        self.n_obj = 2
        self.output_dir = os.path.abspath('../output')
        self.temp_dir = '/tmp/swat_runs'
        self.verbose = True
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

    def __init__(self, cal_config, client=None, **kwargs):
        param = cal_config.cal_param
        n_obj = cal_config.n_obj

        xl = [param[x][0][0] for x in param.keys()]
        xu = [param[x][0][1] for x in param.keys()]

        if client is None:
            ee = True
        else:
            ee = False

        super().__init__(n_var=len(xl),
                         n_obj=n_obj,
                         xl=anp.array(xl),
                         xu=anp.array(xu),
                         evaluation_of=['F'],
                         elementwise_evaluation=ee,
                         **kwargs)

        self.cal_config = cal_config
        self.client = client

    def _evaluate(self, X, out, *args, **kwargs):
        n_obj = self.n_obj
        cal_config = self.cal_config
        client = self.client

        inputs = [[X[k] for k in range(len(X))],
                  [cal_config for k in range(len(X))],
                  [n_obj for k in range(len(X))],
                  [k for k in range(len(X))]]

        if client is not None:
            jobs = client.map(fun, *inputs)
            results = client.gather(jobs)
            f, f_out, param, simulated = zip(*results)

        else:
            f, f_out, param, simulated = fun(X, cal_config, n_obj, 1)

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

def fun(x, cal_config, n_obj, ind):
    swat_model = swat_config.ModelSetup(cal_config.model_file)
    swat_model.agr_treshold = cal_config.agr_treshold
    swat_model.hid_treshold = cal_config.hid_treshold
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
    simulated = run_single_model(swat_model, cal_config.out_var_rch, cal_config.out_var_sub, ind)

    # objective function computation

    q = simulated[0]
    sm = simulated[1]

    start_hid = ['1993-06-01 00:00:00', '1995-01-01 00:00:00', '1997-12-01 00:00:00', '2001-04-01 00:00:00', '2003-01-01 00:00:00', '2009-09-01 00:00:00',\
        '2012-07-01 00:00:00', '2015-05-01 00:00:00']

    end_hid = ['1993-09-01 00:00:00', '1995-04-01 00:00:00', '1998-04-01 00:00:00', '2001-09-01 00:00:00', '2003-03-01 00:00:00', '2010-04-01 00:00:00',\
        '2012-10-01 00:00:00', '2016-01-01 00:00:00']

    start_agr = ['1993-05-01 00:00:00','1994-06-01 00:00:00','1994-11-01 00:00:00','1997-03-01 00:00:00', '2002-01-01 00:00:00','2003-04-01 00:00:00',\
        '2003-12-01 00:00:00', '2006-11-01 00:00:00','2009-11-01 00:00:00','2012-05-01 00:00:00', '2014-07-01 00:00:00','1993-06-01 00:00:00',\
        '1995-01-01 00:00:00', '1997-12-01 00:00:00','2001-04-01 00:00:00','2003-01-01 00:00:00', '2009-09-01 00:00:00','2012-07-01 00:00:00','2015-05-01 00:00:00']

    end_agr = ['1993-09-01 00:00:00','1994-09-01 00:00:00','1995-03-01 00:00:00','1998-03-01 00:00:00', '2002-03-01 00:00:00','2003-06-01 00:00:00',\
       '2004-05-01 00:00:00', '2007-03-01 00:00:00','2010-04-01 00:00:00','2012-10-01 00:00:00', '2016-01-01 00:00:00','1993-09-01 00:00:00',\
       '1995-04-01 00:00:00', '1998-04-01 00:00:00','2001-09-01 00:00:00','2003-03-01 00:00:00', '2010-04-01 00:00:00','2012-10-01 00:00:00','2016-01-01 00:00:00']

    date_start_hid = [datetime.strptime(date, "%Y-%m-%d %H:%M:%S") for date in start_hid]
    date_end_hid = [datetime.strptime(date, "%Y-%m-%d %H:%M:%S") for date in end_hid]

    date_start_agr = [datetime.strptime(date, "%Y-%m-%d %H:%M:%S") for date in start_agr]
    date_end_agr = [datetime.strptime(date, "%Y-%m-%d %H:%M:%S") for date in end_agr]

    factor_list_agr = [0.06,0.06,0.06,0.2,0.06,0.06,0.1,0.06,0.06,0.1,0.2]
    factor_list_hid = [0.1,0.1,0.1,0.1,0.1,0.2,0.1,0.2]

    for factor, datei, datef in zip(factor_list_hid, date_start_hid, date_end_hid):
        f1 = factor*(sum([q_i['severity'].loc[datei:datef].sum() for q_i in q.values()]))

    for factor, datei, datef in zip(factor_list_agr, date_start_agr, date_end_agr):
        f2 = factor*(sum([sm_i['severity'].loc[datei:datef].sum() for sm_i in sm.values()]))

    f = [f1,f2]
    f_out = {}
    f_out['obj_1'] = f1
    f_out['obj_2'] = f2

    return f, f_out, param, simulated

def run_single_model(swat_model, out_var_rch, out_var_sub, ind):
    # assigning ID to swat output folder
    swat_model.new_model_name = 'RUN{:04d}'.format(ind)
    # prepare SWAT run
    swat_model.prepare_swat()
    # execute model
    swat_model.run_swat()
    # get output time series of variables of interest
    ts1 = get_rch_output(swat_model, out_var_rch)
    ts2 = get_sub_output(swat_model, out_var_sub)
    simulated = [ts1, ts2]
    # remove temporal output folders
    swat_model.remove_swat()

    return simulated

def get_sub_output(model, out_var):

    output_dir = model.output_dir + '/' + model.new_model_name + '/output.sub'
    max_col1 = 35
    max_col2 = 35
    ns = 10

    treshold_sm_df = get_sm_drought_treshold(model)
    treshold_group_sub = treshold_sm_df.groupby('sub')

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

        n = 23
        treshold_sub = treshold_group_sub.get_group(sub)
        treshold_sub_se = treshold_sub['treshold']
        treshold_sub_se = pd.concat([treshold_sub_se]*n)

        aux = temp.loc[temp['SUB'] == str(sub), :].copy()
        aux['YEAR'] = years
        aux.index = [dt.datetime(int(year),int(month),1) for year, month in zip(aux['YEAR'], aux['MON'])]

        treshold_sub_se.index = [dt.datetime(int(year),int(month),1) for year, month in zip(aux['YEAR'], aux['MON'])]
        aux['treshold'] = treshold_sub_se
        aux['severity'] = np.where(aux['SWmm']<aux['treshold'], abs(aux['SWmm'] - aux['treshold']), 0)

        output[str(sub)] = {'AREAkm2': aux['AREAkm2'][0], 'Series': aux.loc[:,out_var[0]], 'severity':aux.loc[:,out_var[1]]}

    return output

def get_rch_output(model, out_var):

    output_dir = model.output_dir + '/' + model.new_model_name + '/output.rch'
    max_col1 = 37
    max_col2 = 38
    ns = 12

    treshold_rc_df = get_rc_drought_treshold(model)
    treshold_group_rch = treshold_rc_df.groupby('sub')

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

    temp['sf_mm/s'] = (temp['FLOW_OUTcms'] /(temp['AREAkm2']*1000000))*1000
    temp['sf_mm/d'] = temp['sf_mm/s']*86400

    subbasins = np.sort(np.unique(temp['RCH']).astype(int))

    output = {}
    for sub in subbasins:

        n = 23
        treshold_rch = treshold_group_rch.get_group(sub)
        treshold_rch_se = treshold_rch['treshold']
        treshold_rch_se = pd.concat([treshold_rch_se]*n)

        aux = temp.loc[temp['RCH'] == str(sub), :].copy()
        aux['YEAR'] = years
        aux.index = [dt.datetime(int(year),int(month),1) for year, month in zip(aux['YEAR'], aux['MON'])]

        treshold_rch_se.index = [dt.datetime(int(year),int(month),1) for year, month in zip(aux['YEAR'], aux['MON'])]
        aux['treshold'] = treshold_rch_se
        aux['severity'] = np.where(aux['sf_mm/d']<aux['treshold'], abs(aux['sf_mm/d'] - aux['treshold']), 0)

        output[str(sub)] = {'AREAkm2': aux['AREAkm2'][0], 'Series': aux.loc[:, out_var[0]], 'severity':aux.loc[:,out_var[1]]}

    return output

def get_sm_drought_treshold(model):

    treshold_dir_sm = model.agr_treshold
    col_names = ['gis','month','treshold','sub','drop']
    treshold_sm = open(treshold_dir_sm)

    treshold_sm_df = pd.read_csv(treshold_sm, names = col_names, delimiter= ',', skiprows=1)
    treshold_sm_df['sub'] = treshold_sm_df['sub'].astype(int)

    return treshold_sm_df

def get_rc_drought_treshold(model):

    treshold_dir_rc = model.hid_treshold
    col_names = ['sub','month','treshold','drop1','drop2']
    treshold_rc = open(treshold_dir_rc)

    treshold_rc_df = pd.read_csv(treshold_rc, names = col_names, delimiter= ',', skiprows=1)
    treshold_rc_df['sub'] = treshold_rc_df['sub'].astype(int)

    return treshold_rc_df
