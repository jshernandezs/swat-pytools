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
import subprocess as sp



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
            
            x_g = []
            for array in x:
                x_i = [array[i:i + 7] for i in range(0, len(array), 7)]
                x_g.append(x_i)
            df4 = pd.DataFrame(x_g, columns=df2.columns)
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
        self.simulation_months = os.path.abspath('../resources/csv_files/cienaga_simulation_months_hru.csv.csv')
        self.soils_awc = os.path.abspath('../resources/csv_files/cienaga_awc_soils.csv')
        self.hid_treshold = os.path.abspath('../resources/csv_files/cienaga_hid_dro_treshold.csv')
        self.opt_hrus = os.path.abspath('../resources/csv_files/cienaga_hru_opt.csv')
        self.opt_subs = os.path.abspath('../resources/csv_files/cienaga_sub_opt.csv')
        self.loc_const_1 = os.path.abspath('../resources/csv_files/cienaga_loc_const.csv')

        self.swat_exec_name = 'SWAT_Rev670'
        self.obs_file = ''
        self.out_var_agr = ['sw_mm', 'sw_awc']
        self.out_var_hid = ['q_mm/d', 'q_mm/s']
        self.spatial_units = 933
        self.results_years = 32
        self.pdmm = 7
        self.pdmm_param = {}
        self.n_obj = 2
        self.n_const = 1
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

        n_spatial_units = cal_config.spatial_units
        n_pdmm = cal_config.pdmm
        n_obj = cal_config.n_obj
        n_const = cal_config.n_const

        if client is None:
            ee = True
        else:
            ee = False

        super().__init__(n_var=n_spatial_units*n_pdmm,
                         n_obj=n_obj,
                         n_constr=n_const,
                         xl=0.0,
                         xu=1.0,
                         vtype=int,
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
            f, f_out, g, spatial_unit_opt = zip(*results)

        else:
            f, f_out, g, spatial_unit_opt = fun(X, cal_config, n_obj, 1)

        out["F"] = anp.array(f)
        out["G"] = anp.array(g)

        # update hooker
        lib = Hooker()
        lib.add_sim(f_out, spatial_unit_opt)
        lib.set_output_dir(cal_config.output_dir)

def my_callback(algorithm):
    # write hooked values out
    lib = Hooker()
    lib.print(algorithm)
    lib.clear_hooker()


def fun(x, cal_config, n_obj, ind):

    swat_model = swat_config.ModelSetup(cal_config.model_file)
    swat_model.opt_hrus = cal_config.opt_hrus
    swat_model.opt_subs = cal_config.opt_subs
    swat_model.loc_const_1 = cal_config.loc_const_1

    spatial_unit_ref = dic_hru_sub(swat_model)   
    keysList = list(spatial_unit_ref.keys())
    spatial_unit_opt = {}

    g1 = fun_constraint1(swat_model, x)

    if g1 !=0: 

        for j, key in zip(range(0,len(x),7), enumerate(spatial_unit_ref.keys())):
            spatial_unit_opt[key] = [x[j:j+7]]

        spatial_unit_opt = dict(zip(keysList, list(spatial_unit_opt.values())))

        for key in spatial_unit_opt.keys():
            spatial_unit_opt[key] = [spatial_unit_opt[key][0], spatial_unit_ref[key][0]]

        f1 = float('nan')
        f2 = float('nan')

        f = [f1,f2]

        f_out = {}
        f_out['obj_1'] = f1
        f_out['obj_2'] = f2

    else:    
        swat_model = swat_config.ModelSetup(cal_config.model_file)
        swat_model.opt_hrus = cal_config.opt_hrus
        swat_model.opt_subs = cal_config.opt_subs
        swat_model.simulation_months = cal_config.simulation_months
        swat_model.soils_awc = cal_config.soils_awc
        swat_model.hid_treshold = cal_config.hid_treshold
        swat_model.swat_dir = cal_config.swat_dir
        swat_model.temp_dir = cal_config.temp_run_dir
        swat_model.output_dir = cal_config.temp_dir
        swat_model.verbose = cal_config.verbose

        for j, key in zip(range(0,len(x),7), enumerate(spatial_unit_ref.keys())):
            spatial_unit_opt[key] = [x[j:j+7]]

        spatial_unit_opt = dict(zip(keysList, list(spatial_unit_opt.values())))

        for key in spatial_unit_opt.keys():
            spatial_unit_opt[key] = [spatial_unit_opt[key][0], spatial_unit_ref[key][0]]

        # swat preparation
        swat_model.swat_exec_name = cal_config.swat_exec_name
        swat_model.pdmm_param = cal_config.pdmm_param
        swat_model.spatial_unit = spatial_unit_opt

        # swat execution
        simulated = run_single_model(swat_model, cal_config.out_var_agr, cal_config.out_var_hid, ind)

        # objective function computation
        sm_defict = simulated[0]
        q_deficit = simulated[1]

        f1 = 0
        f2 = 0

        f1 = sum([sm_defict_i['sw_awc'].sum() for sm_defict_i in sm_defict.values()])
        f2 = -1*sum([q_deficit_i['q_mm/d'].sum() for q_deficit_i in q_deficit.values()])

        f = [f1,f2]
        f_out = {}
        f_out['obj_1'] = f1
        f_out['obj_2'] = f2

    return f, f_out, g1, spatial_unit_opt


def run_single_model(swat_model, out_var_agr, out_var_hid, ind):
    # assigning ID to swat output folder
    swat_model.new_model_name = 'RUN{:04d}'.format(ind)
    # prepare SWAT run
    swat_model.prepare_swat()
    # execute model
    swat_model.run_swat()
    # get output time series of variables of interest
    ts1 = get_hru_output(swat_model, out_var_agr)
    ts2 = get_rch_output(swat_model, out_var_hid)
    simulated = [ts1, ts2]
    # remove temporal output folders
    swat_model.remove_swat()

    return simulated

def get_hru_output(model, out_var_agr):

    output_dir = model.output_dir + '/' + model.new_model_name + '/output.hru'

    no_hru = 2500
    results_years = 32

    simulation_months = get_simulation_months(model)
    soils_awc = get_soils_awc(model)

    col_names = range(1,60)
    with open(output_dir) as f:

        df = pd.read_csv(f, names=col_names, sep='\s+',skiprows=9)

        df.reset_index(inplace=True)
        df_f = df.filter(items =['level_0','level_1','level_12'])

    sm_list = []
    for i in range(0,12*no_hru*results_years+(no_hru*results_years),(12*no_hru)+no_hru):
        dfi = df_f.iloc[i:i+12*no_hru]
    
        sm_list.append(dfi)
        sm_df = pd.concat(sm_list, axis=0, ignore_index=True)

    sm_df = sm_df.rename(columns={'level_0':'lu', 'level_1':'hru' ,'level_12':'sw_mm'})

    sm_df['month'] = simulation_months['month']

    sm_df['soil'] = soils_awc['soil']
    sm_df['awc_mm'] = soils_awc['awc_mm']

    land_use_water = ['UWB', 'WATER', 'VAR']
    sm_df = sm_df[sm_df.soil.isin(land_use_water) == False]
    sm_df['awc_mm'] = sm_df['awc_mm'].astype(int)

    wet_months = [8,9,10,11]
    sm_df = sm_df[sm_df.month.isin(wet_months) == False]

    sm_df['sw_awc'] = abs(sm_df['awc_mm']-sm_df['sw_mm'])

    output = {}
    sm_df_hru = sm_df.groupby('hru')

    for name, group in sm_df_hru:

        hru_df = group
        output[str(name)] = {'Series':hru_df.loc[:,out_var_agr[0]], 'sw_awc':hru_df.loc[:,out_var_agr[1]]}
    return output

def get_rch_output(model, out_var_hid):

    output_dir = model.output_dir + '/' + model.new_model_name + '/output.rch'

    no_sub = 205
    results_years = 32
    treshold_rc_df = get_hid_drought_treshold(model)

    col_names = range(1,60)
    with open(output_dir) as f:
        
        df = pd.read_csv(f, names=col_names, sep='\s+',skiprows=9)
        df.reset_index(inplace=True)
        df_f = df.filter([2, 4, 5, 6], axis=1)

    q_list = []
    for i in range(0,12*no_sub*results_years+(no_sub*results_years),(12*no_sub)+no_sub):
        dfi = df_f.iloc[i:i+12*no_sub]
    
        q_list.append(dfi)
        q_df = pd.concat(q_list, axis=0, ignore_index=True)

    q_df = q_df.rename(columns={2: 'sub', 4: 'month', 5: 'area', 6:'q_m3/s'})

    q_df['q_mm/s'] = (q_df['q_m3/s']/(q_df['area']*1000000))*1000
    q_df['q_mm/d'] = q_df['q_mm/s']*86400

    #q_df_aux = q_df.merge(treshold_rc_df, left_on=['sub','month'], right_on=['sub','month'])
    #q_df_aux['severity'] =  abs(q_df_aux['treshold_mm/d']-q_df_aux['q_mm/d'])

    output = {}
    q_df_sub = q_df.groupby('sub')

    for name, group in q_df_sub:

        sub_df = group
        output[str(name)] = {'q_mm/d': sub_df.loc[:, out_var_hid[0]], 'q_mm/s':sub_df.loc[:, out_var_hid[1]]}
    return output

def get_simulation_months(model):

    simulation_months_dir = model.simulation_months
    col_names = ['month']
    simulation_month = open(simulation_months_dir)

    simulation_months_df = pd.read_csv(simulation_month, names = col_names, delimiter= ',', skiprows=1)
    return simulation_months_df 

def get_soils_awc(model):

    soils_awc_dir = model.soils_awc
    col_names = ['soil','awc_mm']
    soil_awc = open(soils_awc_dir)

    soils_awc_df = pd.read_csv(soil_awc, names = col_names, delimiter= ',', skiprows=1)
    return soils_awc_df

def get_hid_drought_treshold(model):

    treshold_hid_dir = model.hid_treshold
    col_names = ['sub','month','treshold_mm/d']
    treshold_hid = open(treshold_hid_dir)

    treshold_rc_df = pd.read_csv(treshold_hid, names = col_names, delimiter= ',', skiprows=1)
    treshold_rc_df['sub'] = treshold_rc_df['sub'].astype(int)
    return treshold_rc_df

def dic_hru_sub(model):

    opt_hrus_dir = model.opt_hrus
    col_names = ['hru', 'soil']
    opt_hrus = open(opt_hrus_dir)

    opt_hrus_df = pd.read_csv(opt_hrus, names = col_names, delimiter= ',', skiprows=1, converters={'hru': str})
    hrus, soil  = list(opt_hrus_df['hru']), list(opt_hrus_df['soil'])

    hrus_dict = {z[0]: list(z[1:]) for z in zip(hrus, soil)}

    opt_sub_dir = model.opt_subs
    col_names = ['sub', 'ext']
    opt_sub = open(opt_sub_dir)

    opt_subs_df = pd.read_csv(opt_sub, names = col_names, delimiter= ',', skiprows=1)
    subs, value  = list(opt_subs_df['sub']), list(opt_subs_df['ext'])

    subs_dict = {z[0]: list(z[1:]) for z in zip(subs, value)}

    spatial_unit_ref = {**hrus_dict, **subs_dict}
    return spatial_unit_ref
   
def fun_constraint1(model, x):

    const_dir = model.loc_const_1
    const_file = open(const_dir)
    loc_const = anp.genfromtxt(const_file, delimiter=';').flatten()

    g1 = anp.sum(x*loc_const, axis=0)
    return g1


  
