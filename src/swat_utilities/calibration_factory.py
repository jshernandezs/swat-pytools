#!/usr/bin/env python3
""" Factory of problems for SWAT model calibration
"""
import os
import autograd.numpy as anp
from pymoo.model.problem import Problem
from pymoo.model.individual import Individual
from pymoo.model.population import Population
from pymoo.model.algorithm import filter_optimum
from swat_utilities import swat_config
from swat_utilities.performance import Metric, read_obs_file
from swat_utilities.hit import HydrologicIndex as Hi
from multiprocessing import Lock
import pandas as pd
import numpy as np
import subprocess as sp


class Hooker:
    """ Stores simulations for each generation/iteration
    """
    simulations = []
    params = []
    obj_f = []
    output_dir = ''

    def add_sim(self, f, param, sim):

        new_param = {key: param[key][0] for key in param.keys()}

        self.obj_f.append(f)
        self.simulations.append(sim)
        self.params.append(new_param)

    def print(self, algorithm):
        n_gen = algorithm.n_gen
        if not self.output_dir == '':
            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)

            df1 = pd.concat(self.simulations, axis=1)
            df2 = pd.DataFrame(self.params)
            df3 = pd.DataFrame(self.obj_f)

            filename1 = 'simulations_gen{:04d}.csv'.format(n_gen)
            filename2 = 'parameters_gen{:04d}.csv'.format(n_gen)
            filename3 = 'objs_gen{:04d}.csv'.format(n_gen)

            df1.to_csv(self.output_dir + '/' + filename1, header=False, index=True)
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
        Hooker.simulations = []
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
        self.out_file = 'watout.dat'
        self.out_var = 'FLOWm^3/s'
        self.cal_param = {}
        self.obj_f = {'NSE': 'none'}
        self.constrs = {}
        self.output_dir = os.path.abspath('../output')
        self.temp_dir = '/tmp/swat_runs'
        self.verbose = True
        self.lock = 0
        self.runid = swat_config.RunID()
        self.temp_run_dir = '/tmp/output_swat'
        self.rds_file = ''
        self.opt_ncv = False
        self.ind_dict = {}

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

        obj_f = cal_config.obj_f
        aux = [[item] if type(item) in [str, int, float] else item for item in obj_f.values()]
        n_obj = len([item for sublist in aux for item in sublist])

        constrs = cal_config.constrs
        n_constr = len(constrs.values())

        xl = [param[x][0][0] for x in param.keys()]
        xu = [param[x][0][1] for x in param.keys()]

        super().__init__(n_var=len(xl),
                         n_obj=n_obj,
                         n_constr=n_constr,
                         xl=anp.array(xl),
                         xu=anp.array(xu),
                         evaluation_of=['F', 'G', 'CV', 'NCV'],
                         elementwise_evaluation=True,
                         **kwargs)

        cal_config.lock = Lock()

        self.cal_config = cal_config

    def _evaluate(self, x, out, *args, **kwargs):
        n_obj = self.n_obj
        n_constr = self.n_constr
        cal_config = self.cal_config

        f, g, f_out, param, simulated = fun(x, cal_config, n_obj, n_constr)
        out["F"] = anp.array(f)
        out["G"] = anp.array(g)

        # update hooker
        lib = Hooker()
        lib.add_sim(f_out, param, simulated)
        lib.set_output_dir(cal_config.output_dir)


def my_callback(algorithm):
    # write hooked values out
    lib = Hooker()
    lib.print(algorithm)
    lib.clear_hooker()
    # clear list of ids
    ids = swat_config.RunID()
    ids.clear_id()


def fun(x, cal_config, n_obj, n_constr):
    lock = cal_config.lock
    runid = cal_config.runid
    obj_f = cal_config.obj_f
    constrs = cal_config.constrs
    opt_ncv = cal_config.opt_ncv
    ind_dict = cal_config.ind_dict
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
    simulated = swat_config.run_single_model(swat_model, cal_config.out_file, cal_config.out_var, lock, runid)

    # objective function computation
    obs_file = cal_config.obs_file
    simulation = Metric(obs_file, simulated)

    transforms = [[item] if isinstance(item, str) else item for item in obj_f.values()]

    f = anp.empty([1, n_obj])
    f_out = {}
    f_eval = []
    idx = 0

    objs = list(obj_f.keys())
    perf_list = [obj for obj in objs if obj in ['NSE', 'KGE', 'KGEp', 'IoA', 'R2', 'koR2', 'ksR2', 'RMSE', 'R4MS4E',
                                                'RSR', 'PBIAS']]
    t_new = [t for t, obj in zip(transforms, objs) if obj in ['NSE', 'KGE', 'KGEp', 'IoA', 'R2', 'koR2', 'ksR2', 'RMSE',
                                                              'R4MS4E', 'RSR', 'PBIAS']]
    for i, obj in enumerate(perf_list):
        for transf in t_new[i]:
            if obj in ['NSE', 'KGE', 'KGEp', 'IoA', 'R2', 'koR2', 'ksR2']:
                f_eval = 1 - getattr(simulation, 'get_{:s}'.format(obj.lower()))(transf.lower())
            elif obj in ['RMSE', 'R4MS4E', 'RSR']:
                f_eval = getattr(simulation, 'get_{:s}'.format(obj.lower()))(transf.lower())
            elif obj == 'PBIAS':
                f_eval = abs(simulation.get_pbias())

            name = '{:s}_{:s}'.format(obj, transf)
            f_out[name] = f_eval
            f[0, idx] = f_eval
            idx += 1

    hi_list = [obj for obj in objs if obj not in ['NSE', 'KGE', 'KGEp', 'IoA', 'R2', 'koR2', 'ksR2', 'RMSE',
                                                  'R4MS4E', 'RSR', 'PBIAS']]
    threshs = [t for t, obj in zip(transforms, objs) if obj not in ['NSE', 'KGE', 'KGEp', 'IoA', 'R2', 'koR2', 'ksR2', 
                                                                    'RMSE', 'R4MS4E', 'RSR', 'PBIAS']]
    if len(hi_list) > 0:
        eps = np.finfo(float).eps
        observed = read_obs_file(obs_file)
        obs_hit = Hi(observed)
        sim_hit = Hi(simulated, obs_file=obs_file)

        for i, hi_axis in enumerate(hi_list):
            ref = threshs[i]
            inds = ind_dict[hi_axis]
            index_obs_dict = obs_hit.get_list_indices(inds)
            index_sim_dict = sim_hit.get_list_indices(inds)

            error = []
            for ind in inds:
                index_obs = index_obs_dict[ind.upper()]
                index_sim = index_sim_dict[ind.upper()]

                error_i = abs(index_obs - index_sim) / abs(index_obs + eps) * 100 - ref
                error.append(error_i)

            error = np.array(error)
            f_eval = np.nanmean(error[error >= 0])
            name = '{:s}_mean_rel_abs_error'.format(hi_axis)
            f_out[name] = f_eval
            f[0, idx] = f_eval
            idx += 1

    # constraints computation

    g = []

    if n_constr > 0:

        g_out = {}
        idx = 0

        observed = read_obs_file(obs_file)
        obs_hit = Hi(observed)

        perf_list = ['NSE', 'KGE', 'KGEp', 'IoA', 'R2', 'koR2', 'ksR2', 'RMSE', 'R4MS4E', 'RSR', 'PBIAS']
        perf_constr = [constr for constr in constrs.keys() if constr in perf_list]
        hi_constr = [constr for constr in constrs.keys() if (constr in obs_hit.list_all) or
                     (constr in obs_hit.list_mag_seven)]

        nc = len(perf_constr)
        if len(hi_constr) > 0:
            nc = nc + 1

        g = anp.empty([1, nc])
        g_eval = -9999
        name = ''

        for constr in perf_constr:

            if constr in ['NSE', 'KGE', 'KGEp', 'IoA', 'R2', 'koR2', 'ksR2']:
                transf = constrs[constr][0]
                thresh = constrs[constr][1]
                g_eval = (1 - getattr(simulation, 'get_{:s}'.format(constr.lower()))(transf.lower())) / thresh - 1
                name = 'constr_{:s}_{:s}_{:.1f}'.format(constr, transf, thresh)
                
            elif constr in ['RMSE', 'R4MS4E', 'RSR']:
                transf = constrs[constr][0]
                thresh = constrs[constr][1]
                g_eval = getattr(simulation, 'get_{:s}'.format(constr.lower()))(transf.lower())
                name = 'constr_{:s}_{:s}_{:.1f}'.format(constr, transf, thresh)

            elif constr == 'PBIAS':
                thresh = constrs[constr]
                g_eval = abs(simulation.get_pbias()) / thresh - 1
                name = 'constr_{:s}_{:.1f}'.format(constr, thresh)

            g_out[name] = g_eval
            g[0, idx] = g_eval
            idx += 1

        if len(hi_constr) > 0:
            eps = np.finfo(float).eps
            sim_hit = Hi(simulated, obs_file=obs_file)
            index_obs_dict = obs_hit.get_list_indices(hi_constr)
            index_sim_dict = sim_hit.get_list_indices(hi_constr)

            cv = []
            w = []
            for constr in hi_constr:
                index_obs = index_obs_dict[constr.upper()]
                index_sim = index_sim_dict[constr.upper()]

                if type(constrs[constr]) == list:
                    thresh = constrs[constr][0]
                else:
                    thresh = constrs[constr]

                wi = 1/len(hi_constr)

                if (type(constrs[constr]) in [list]) and (len(constrs[constr]) > 1):
                    wi = constrs[constr][1]

                g_eval = abs(index_obs - index_sim) / abs(index_obs + eps) * 100 / thresh - 1
                name = 'constr_{:s}_{:.1f}'.format(constr, thresh)

                g_out[name] = g_eval
                cv.append(g_eval)
                w.append(wi)

            cv = np.array(cv)
            w = np.array(w)

            w = w[cv > 0]
            cv = cv[cv > 0]

            cv_eq = sum([wi * cvi for wi, cvi in zip(w, cv)])

            if opt_ncv:
                # nw = [round(wi*1000, 0) if wi != 1 else wi for wi in w]
                # ncv = sum(nw)
                ncv = cv.size
                cv_eq = ncv + cv_eq

            g[0, idx] = cv_eq

        f_out.update(g_out)

    return f, g, f_out, param, simulated
