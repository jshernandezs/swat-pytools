#!/usr/bin/env python3
"""
Classes for processing SWAT calibration results
"""
import pandas as pd
import numpy as np
import datetime as dt
import os
import matplotlib.pyplot as plt
from swat_utilities.performance import read_obs_file


class Optimal:
    def __init__(self, cal_dir):
        dir_list = os.listdir(cal_dir)      
        dir_list.sort()   

        self.cal_dir = cal_dir
        self.objs_list = [f for f in dir_list if f.startswith('opt_objs_')]
        self.parameter_list = [f for f in dir_list if f.startswith('opt_parameters_')]
        self.sim_list = [f for f in dir_list if f.startswith('simulations_')]
        self.param_ref_list = [f for f in dir_list if f.startswith('parameters_')]
        self.obj_ref_list = [f for f in dir_list if f.startswith('objs_')]
        self.ngen = []
        self.param_names = []
        self.best_of = pd.DataFrame()
        self.best_parameters = pd.DataFrame()
        self.best_simulation = pd.DataFrame()
        self.address = []
        
    def get_obj(self, ngen=None):
        
        cal_dir = self.cal_dir
        if ngen is None:
            ngen = len(self.objs_list)
        self.ngen = ngen
        
        objs_file = self.objs_list[ngen - 1]
        self.best_of = pd.read_csv(cal_dir + '/' + objs_file, header=0)
        
        return self.best_of
    
    def get_params(self, ngen=None):
        
        cal_dir = self.cal_dir
        if ngen is None:
            ngen = len(self.objs_list)
        self.ngen = ngen
        
        parameter_file = self.parameter_list[ngen - 1]
        self.best_parameters = pd.read_csv(cal_dir + '/' + parameter_file, header=0)
        
        return self.best_parameters

    def get_best(self, ngen=None):
        
        cal_dir = self.cal_dir
        if ngen is None:
            ngen = len(self.objs_list)
        self.ngen = ngen

        objs_file = self.objs_list[ngen - 1]
        parameter_file = self.parameter_list[ngen - 1]
        sim_list = self.sim_list[:ngen]
        param_ref_list = self.param_ref_list[:ngen]
        
        print('Obtaining time-series associated to Pareto-optimal solutions...')
        address = []
        aux = pd.read_csv(cal_dir + '/' + parameter_file).to_numpy()
        index = np.arange(0, aux.shape[0])
        c = 1
        while 1:
            ref_params = pd.read_csv(cal_dir + '/' + param_ref_list[-c]).to_numpy()
            with open(cal_dir + '/'  + sim_list[-c], 'rb') as f:
                flow = pd.read_csv(f, header=None, index_col=0)
                
            dates = [dt.datetime.strptime(t, '%Y-%m-%d') for t in flow.index]
            flow.index = dates
    
            gen = int(param_ref_list[-c][-8:-4])
            ind = []
            for j, item in enumerate(aux):
                for i, ref in enumerate(ref_params):
                    if (np.round(item, 7) == np.round(ref, 7)).all():
                        ts = flow.iloc[:, i]
                        address.append([index[j], gen, i, ts])
                        ind.append(j)
            aux = np.delete(aux, ind, 0)
            index = np.delete(index, ind)
            if aux.size == 0:
                break
            else:
                c += 1
        
        address = np.array(address, dtype=object)
        address = address[address[:,0].argsort()]
        
        best_simulation = address[:,-1].tolist()
        parameters = pd.read_csv(cal_dir + '/' + parameter_file, header=0)
                
        self.param_names = parameters.columns
        self.best_of = pd.read_csv(cal_dir + '/' + objs_file, header=0)
        self.best_parameters = parameters
        self.address = list(address[:,:-1])
        try:
            self.best_simulation = pd.concat(best_simulation, axis=1)
        except ValueError:
            self.best_simulation = best_simulation

        return self.best_of, self.best_parameters, self.best_simulation

    def write_best(self, output_dir, obs_file=None):

        if (output_dir is not None) & (len(self.best_of) > 0):
            try:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)

                best_of = self.best_of
                best_sim = self.best_simulation
                best_parameters = self.best_parameters

                if len(obs_file) > 0:
                    try:
                        obs = read_obs_file(obs_file)
                        best_sim = best_sim[obs.index[0]:obs.index[-1]]
                    except ValueError:
                        print('Something went wrong with the observed file')

                best_of.to_csv(output_dir + '/' + 'best_functions_gen{:04d}.csv'.format(self.ngen),
                               header=True, index=False)
                best_parameters.to_csv(output_dir + '/' + 'best_parameters_gen{:04d}.csv'.format(self.ngen),
                                       header=True, index=False)
                best_sim.to_csv(output_dir + '/' + 'best_simulation_gen{:04d}.csv'.format(self.ngen),
                                header=False, index=True)

            except OSError:
                print('Something went wrong with the output directory')
        else:
            print('You must run get_best() method first')

    def plot_ts(self, obs_file=None, opt_log=False):

        sim = self.best_simulation
        try:   
            fig, ax = plt.subplots()

            if obs_file is not None:
                try:
                    obs = read_obs_file(obs_file)
                    sim = sim[obs.index[0]:obs.index[-1]]
                    ax.plot(obs.index, obs.values, color='b', label='observed', linewidth=0.7, zorder=2)
                    ax.plot(sim.index, sim.values, color='grey', label='simulated', zorder=1)
                    ax.legend(ax.lines[:2], ['observed', 'simulated'])
                except ValueError:
                    print('Something went wrong with the observed file')
            else:
                ax.plot(sim.index, sim.values, color='grey', label='simulated')
                # ax.legend(ax.lines[0], ['simulated'])
            
            plt.ylabel('Flow ($m^3 s^{-1}$)')
                
            if opt_log:
                plt.yscale('log')
                plt.grid(True, which="major", linewidth=1)
                plt.grid(True, which="minor", linestyle='--')

            plt.show()

        except RuntimeError:
            print('You must run get_best() method first')
    
    def plot_fdc(self, obs_file=None, opt_log=True):
        
        sim = self.best_simulation        
        try:   
            fig, ax = plt.subplots()

            if obs_file is not None:
                try:
                    obs = read_obs_file(obs_file)
                    sim = sim[obs.index[0]:obs.index[-1]]
                    
                    q_sim = fdc(sim)
                    q_obs = fdc(obs)
                    
                    ax.plot(q_obs, color='b', label='observed', zorder=2)
                    ax.plot(q_sim, color='grey', label='simulated', zorder=1)
                    ax.legend(ax.lines[:2], ['observed', 'simulated'])
                except ValueError:
                    print('Something went wrong with the observed file')
            else:
                q_sim = fdc(sim)
                ax.plot(q_sim, color='grey', label='simulated')
                # ax.legend(ax.lines[0], ['simulated'])
            
            plt.ylabel('Flow ($m^3 s^{-1}$)')
            plt.xlabel('Exceedance Probability (%)')
            
            plt.grid(True, which="major", linewidth=1)
            plt.grid(True, which="minor", linestyle='--')
                
            if opt_log:
                plt.yscale('log')                

            plt.show()

        except RuntimeError:
            print('You must run get_best() method first')
        

def fdc(values):
    values = np.array(values)
    ncol = values.shape[1]
    values_sorted = []
    for i in range(0, ncol):
        # create numpy array for discharge values
        val = values[:, i]
        # sort values
        values_sorted.append(np.sort(val))
        
    # probabilities
    n = values.shape[0]
    p = 1 - np.linspace(1, n, n) / (1 + n)
    
    values_sorted = pd.DataFrame(values_sorted).transpose()
    values_sorted.index = p*100
        
    return values_sorted

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", title="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.
    Source: https://matplotlib.org/
    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="3%", pad=0.05)
    # cbar = ax.figure.colorbar(im, ax=ax, cax=cax, **cbar_kw)
    # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels, fontsize=10)
    ax.set_yticklabels(row_labels, fontsize=10)

    # Let the horizontal axes labeling appear on bottom.
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=40, ha="right",
             va="top", rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_xlabel(title, fontsize=18)

    return im#, cbar