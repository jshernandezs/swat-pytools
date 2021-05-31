#!/usr/bin/env python3
"""
Class containing set of metrics for model evaluation
"""
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class Metric:
    def __init__(self, obs_file, simulated):

        observed = read_obs_file(obs_file)
        series = pd.merge(observed, simulated, how='inner', left_index=True, right_index=True)
        series.columns = ['observed', 'simulated']

        self.series = series
        self.obs = observed
        self.sim = simulated

        self.NSE = []
        self.RMSE = []
        self.KGE = []
        self.KGEp = []
        self.RSR = []
        self.IoA = []
        self.R2 = []
        self.koR2 = []
        self.ksR2 = []
        self.R4MS4E = []
        self.PBIAS = []
        self.ks = []
        self.ko = []
        self.alpha = []
        self.beta = []
        self.gamma = []
        self.r = []

        self.NSE_transform = ''
        self.RMSE_transform = ''
        self.KGE_transform = ''
        self.KGEp_transform = ''
        self.RSR_transform = ''
        self.IoA_transform = ''
        self.R2_transform = ''
        self.koR2_transform = ''
        self.ksR2_transform = ''
        self.R4MS4E_transform = ''
        self.PBIAS_transform = ''
        self.comps_transform = ''

    def get_all(self, transform='none'):

        self.get_nse(transform)
        self.get_rmse(transform)
        self.get_r4ms4e()
        self.get_kge(transform)
        self.get_kgep(transform)
        self.get_rsr(transform)
        self.get_pbias()
        self.get_r2(transform)
        self.get_ioa(transform)
        self.get_kor2(transform)
        self.get_ksr2(transform)
        self.get_kge_components(transform)

        return (self.NSE, self.RMSE, self.KGE, self.KGEp, self.RSR, self.IoA, self.R2, self.koR2, self.ksR2, 
                self.R4MS4E, self.PBIAS)

    def get_nse(self, transform='none'):
        obs0 = np.array(self.series['observed'])
        sim0 = np.array(self.series['simulated'])

        obs, sim = trans(obs0, sim0, opt=transform)

        if transform == 'rel':
            den1 = obs
            den2 = obs.mean()
        else:
            den1 = 1
            den2 = 1

        res = (sim - obs) / den1
        obs_mean = obs.mean()
        res_obs = (obs - obs_mean) / den2

        self.NSE = 1 - res.dot(res.T) / res_obs.dot(res_obs.T)
        self.NSE_transform = transform

        return self.NSE

    def get_rsr(self, transform='none'):
        obs0 = np.array(self.series['observed'])
        sim0 = np.array(self.series['simulated'])

        obs, sim = trans(obs0, sim0, opt=transform)

        if transform == 'rel':
            den1 = obs
            den2 = obs.mean()
        else:
            den1 = 1
            den2 = 1

        res = (sim - obs) / den1
        obs_mean = obs.mean()
        res_obs = (obs - obs_mean) / den2

        nse = 1 - res.dot(res.T) / res_obs.dot(res_obs.T)
        self.RSR = np.sqrt(1 - nse)
        self.RSR_transform = transform

        return self.RSR

    def get_pbias(self):
        obs = np.array(self.series['observed'])
        sim = np.array(self.series['simulated'])

        res = sim - obs
        self.PBIAS = res.sum() / obs.sum() * 100

        return self.PBIAS

    def get_rmse(self, transform='none'):
        obs0 = np.array(self.series['observed'])
        sim0 = np.array(self.series['simulated'])

        obs, sim = trans(obs0, sim0, opt=transform)

        if transform == 'rel':
            den = obs
        else:
            den = 1

        res = (sim - obs) / den
        self.RMSE = np.sqrt(res.dot(res.T) / len(res))
        self.RMSE_transform = transform

        return self.RMSE

    def get_kge(self, transform='none', sr=1, sa=1, sb=1):
        obs0 = np.array(self.series['observed'])
        sim0 = np.array(self.series['simulated'])

        obs, sim = trans(obs0, sim0, opt=transform)

        if transform == 'rel':
            den1 = obs.mean()
            den2 = sim.mean()
        else:
            den1 = 1
            den2 = 1

        obs_mean = obs.mean()
        sim_mean = sim.mean()

        res_obs = (obs - obs_mean) / den1
        res_sim = (sim - sim_mean) / den2
        n = len(obs)

        sigma_obs = np.sqrt(res_obs.dot(res_obs.T) / n)
        sigma_sim = np.sqrt(res_sim.dot(res_sim.T) / n)
        cov = res_sim.dot(res_obs.T) / n

        r = cov / (sigma_obs * sigma_sim)
        alpha = sigma_sim / sigma_obs
        beta = sim_mean / obs_mean

        ed = np.sqrt(sr * (r - 1) ** 2 + sa * (alpha - 1) ** 2 + sb * (beta - 1) ** 2)
        self.KGE = 1 - ed
        self.KGE_transform = transform

        return self.KGE
    
    def get_kgep(self, transform='none', sr=1, sa=1, sb=1):
        obs0 = np.array(self.series['observed'])
        sim0 = np.array(self.series['simulated'])

        obs, sim = trans(obs0, sim0, opt=transform)

        if transform == 'rel':
            den1 = obs.mean()
            den2 = sim.mean()
        else:
            den1 = 1
            den2 = 1

        obs_mean = obs.mean()
        sim_mean = sim.mean()

        res_obs = (obs - obs_mean) / den1
        res_sim = (sim - sim_mean) / den2
        n = len(obs)

        sigma_obs = np.sqrt(res_obs.dot(res_obs.T) / n)
        sigma_sim = np.sqrt(res_sim.dot(res_sim.T) / n)
        cov = res_sim.dot(res_obs.T) / n

        r = cov / (sigma_obs * sigma_sim)
        gamma = (sigma_sim * obs_mean) / (sigma_obs * sim_mean)
        beta = sim_mean / obs_mean

        ed = np.sqrt(sr * (r - 1) ** 2 + sa * (gamma - 1) ** 2 + sb * (beta - 1) ** 2)
        self.KGEp = 1 - ed
        self.KGEp_transform = transform

        return self.KGEp
    
    def get_kge_components(self, transform='none'):
        
        obs0 = np.array(self.series['observed'])
        sim0 = np.array(self.series['simulated'])

        obs, sim = trans(obs0, sim0, opt=transform)

        if transform == 'rel':
            den1 = obs.mean()
            den2 = sim.mean()
        else:
            den1 = 1
            den2 = 1

        obs_mean = obs.mean()
        sim_mean = sim.mean()

        res_obs = (obs - obs_mean) / den1
        res_sim = (sim - sim_mean) / den2
        n = len(obs)

        sigma_obs = np.sqrt(res_obs.dot(res_obs.T) / n)
        sigma_sim = np.sqrt(res_sim.dot(res_sim.T) / n)
        cov = res_sim.dot(res_obs.T) / n

        r = cov / (sigma_obs * sigma_sim)
        alpha = sigma_sim / sigma_obs
        gamma = (sigma_sim * obs_mean) / (sigma_obs * sim_mean)
        beta = sim_mean / obs_mean
        
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.r = r
        self.comps_transform = transform        

    def get_r2(self, transform='none'):
        obs0 = np.array(self.series['observed'])
        sim0 = np.array(self.series['simulated'])

        obs, sim = trans(obs0, sim0, opt=transform)

        if transform == 'rel':
            den1 = obs.mean()
            den2 = sim.mean()
        else:
            den1 = 1
            den2 = 1

        obs_mean = obs.mean()
        sim_mean = sim.mean()

        res_obs = (obs - obs_mean) / den1
        res_sim = (sim - sim_mean) / den2
        n = len(obs)

        sigma_obs = np.sqrt(res_obs.dot(res_obs.T) / n)
        sigma_sim = np.sqrt(res_sim.dot(res_sim.T) / n)
        cov = res_sim.dot(res_obs.T) / n

        r = cov / (sigma_obs * sigma_sim)
        self.R2 = r ** 2
        self.R2_transform = transform

        return self.R2

    def get_kor2(self, transform='none'):
        obs0 = np.array(self.series['observed'])
        sim0 = np.array(self.series['simulated'])

        obs, sim = trans(obs0, sim0, opt=transform)

        if transform == 'rel':
            den1 = obs.mean()
            den2 = sim.mean()
        else:
            den1 = 1
            den2 = 1

        obs_mean = obs.mean()
        sim_mean = sim.mean()

        res_obs = (obs - obs_mean) / den1
        res_sim = (sim - sim_mean) / den2
        n = len(obs)

        sigma_obs = np.sqrt(res_obs.dot(res_obs.T) / n)
        sigma_sim = np.sqrt(res_sim.dot(res_sim.T) / n)
        cov = res_sim.dot(res_obs.T) / n

        r = cov / (sigma_obs * sigma_sim)

        y = sim.reshape(n, 1) / den2  # simulations in the Y-axis
        x = np.concatenate((np.ones((n, 1)), obs.reshape(n, 1) / den1), axis=1)

        try:
            beta = np.linalg.solve(x.T.dot(x), x.T.dot(y))
        except np.LinAlgError:
            beta = [float('inf'), 0]  # assuming singularity problem

        b = beta[1]
        if abs(b) < 1:
            b_r2 = abs(b) * r ** 2
        else:
            b_r2 = (r ** 2) / abs(b)

        self.koR2 = float(b_r2)
        self.ko = beta
        self.koR2_transform = transform

        return self.koR2

    def get_ksr2(self, transform='none'):
        obs0 = np.array(self.series['observed'])
        sim0 = np.array(self.series['simulated'])

        obs, sim = trans(obs0, sim0, opt=transform)

        if transform == 'rel':
            den1 = obs.mean()
            den2 = sim.mean()
        else:
            den1 = 1
            den2 = 1

        obs_mean = obs.mean()
        sim_mean = sim.mean()

        res_obs = (obs - obs_mean) / den1
        res_sim = (sim - sim_mean) / den2
        n = len(obs)

        sigma_obs = np.sqrt(res_obs.dot(res_obs.T) / n)
        sigma_sim = np.sqrt(res_sim.dot(res_sim.T) / n)
        cov = res_sim.dot(res_obs.T) / n

        r = cov / (sigma_obs * sigma_sim)

        y = obs.reshape(n, 1) / den1  # observations in the Y-axis
        x = np.concatenate((np.ones((n, 1)), sim.reshape(n, 1) / den2), axis=1)

        try:
            beta = np.linalg.solve(x.T.dot(x), x.T.dot(y))
        except np.LinAlgError:
            beta = [float('inf'), 0]  # assuming singularity problem

        b = beta[1]
        if abs(b) < 1:
            b_r2 = abs(b) * r ** 2
        else:
            b_r2 = (r ** 2) / abs(b)

        self.ksR2 = float(b_r2)
        self.ks = beta
        self.ksR2_transform = transform

        return self.ksR2

    def get_r4ms4e(self, transform='none'):
        obs0 = np.array(self.series['observed'])
        sim0 = np.array(self.series['simulated'])

        obs, sim = trans(obs0, sim0, opt=transform)

        if transform == 'rel':
            den = obs
        else:
            den = 1

        res = (sim - obs) / den
        num = np.sum(np.power(res, 4))

        self.R4MS4E = np.power(num / len(res), 0.25)
        return self.R4MS4E

    def get_ioa(self, transform='none'):
        obs0 = np.array(self.series['observed'])
        sim0 = np.array(self.series['simulated'])

        obs, sim = trans(obs0, sim0, opt=transform)

        if transform == 'rel':
            den1 = obs
            den2 = obs.mean()
        else:
            den1 = 1
            den2 = 1

        res = (sim - obs) / den1
        obs_mean = obs.mean()
        res_obs = (obs - obs_mean) / den2
        res_sim = (sim - obs_mean) / den2
        res_mix = np.abs(res_obs) + np.abs(res_sim)

        self.IoA = 1 - res.dot(res.T) / res_mix.dot(res_mix.T)
        self.IoA_transform = transform

        return self.IoA

    def plot(self):

        obs = self.obs
        sim = self.sim

        date_obs = obs.index
        sim = sim[date_obs[0]:date_obs[-1]]

        styles = ['r--', 'k-']
        fig, ax = plt.subplots()
        obs.plot(style=styles[0], ax=ax)
        sim.plot(style=styles[1], ax=ax)
        plt.gca().legend(('observed', 'simulated'))
        plt.show()

    def write_obs(self, filename):
        obs = self.obs
        obs.to_csv(filename, header=False)

    def write_sim(self, filename):
        sim = self.sim
        sim.to_csv(filename, header=False)


def trans(obs, sim, opt='none'):
    eps = 1e-6

    if opt == 'log':
        obs2 = np.log(obs + eps)
        sim2 = np.log(sim + eps)
    elif opt == 'inverse':
        obs2 = np.reciprocal(obs + eps)
        sim2 = np.reciprocal(sim + eps)
    elif opt == 'sqrt':
        obs2 = np.sqrt(obs)
        sim2 = np.sqrt(sim)
    else:
        obs2 = obs + eps
        sim2 = sim + eps

    return obs2, sim2


def read_obs_file(obs_file):
    with open(obs_file) as f:
        lines = f.readlines()
        aux = [y.split(',') for y in [x.strip() for x in lines]]
        years = [int(x[0]) for x in aux]
        days = [int(x[1]) for x in aux]
        values = [float(x[2]) for x in aux]

    dates = [dt.datetime(year, 1, 1) + dt.timedelta(day - 1) for year, day in zip(years, days)]
    observed = pd.DataFrame(values, index=dates, columns=['Values'])
    observed = observed.dropna()
    return observed
