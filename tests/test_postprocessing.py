#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import hvwfg
import pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
import matplotlib.patches as patches
from swat_utilities.performance import read_obs_file
from swat_utilities.post_processing import Optimal, fdc
from swat_utilities.performance import Metric
from swat_utilities.hit import HydrologicIndex as Hi
from pymoo.factory import get_decision_making


def main(obs_file, obs_file_val, list_ind, groups, ind_dict, perf_dict, cal_dirs, max_gen, val_file, opt = [True, True, True], 
         output_dir = '../figures'):
       
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Get hypervolume and number of  Pareto solutions over generations

    hv_list, npoints_list, best_list, feasible_gen = get_hv_list(cal_dirs, max_gen)
  
    # Get initial population
    
    swat_optimal = Optimal(cal_dirs[0])
    swat_optimal.objs_list = swat_optimal.obj_ref_list
    swat_optimal.parameter_list = swat_optimal.param_ref_list
    initial = swat_optimal.get_best(ngen=1)
    best_list.append(initial)
    
    # Generate performance plot
    
    if opt[0]:    
        plot_performance(hv_list, npoints_list, best_list, feasible_gen, output_dir)
    
    # Prepare csv files for plotting boxplots in R
    
    with open(val_file, 'rb') as file:
        val_list = pickle.load(file)
    
    if opt[1]:
        all_errors_cal = get_ERHI_performance(obs_file, best_list, list_ind, groups, output_dir)
        all_errors_val = get_ERHI_performance(obs_file_val, val_list, list_ind, groups, output_dir, opt_val=True)
        
    # get FDC and water balance signatures
    fdc_wb_sig, metrics = get_stats(obs_file, best_list, ind_dict, perf_dict)
    fdc_wb_sig_val, metrics_val = get_stats(obs_file_val, val_list, ind_dict, perf_dict, opt_val=True)
    
    # plot PCP comparing objective functions
    # cols = list(range(0,12))
    # limits = [[0.5,1],[-1.5,1],[0,1],[0,1],[0,1],[0,1],[0,30],[0,45],[0,30],[0,120],[0,30],[0,30]]
    # ynames = ['$KGE$', '$KGE_{inv}$', '$R^2$', '$R^2_{sqrt}$', '$R^2_{inv}$', '$IoA_{rel}$',
    #           '$IHA_{1}$', '$IHA_{2}$', '$IHA_{3}$', '$IHA_{4}$', '$IHA_{5}$', 'MAG']    
    # plot_PCP(metrics, cols, limits, ynames, alpha=0.5)
    
    # plot FDCs for best tradeoffs
       
    if opt[2]:
        ind_tradeoff = plot_tradeoffs(obs_file, obs_file_val, best_list, val_list, 
                                      fdc_wb_sig, fdc_wb_sig_val, output_dir)
    
    # get summary of Pareto-optimal performance
    po_stats = get_perf_stats(metrics, ind_tradeoff)
    po_stats_val = get_perf_stats(metrics_val, ind_tradeoff)
    
    return fdc_wb_sig, fdc_wb_sig_val, po_stats, po_stats_val, all_errors_cal, all_errors_val

def get_perf_stats(metrics, ind_tradeoff):
    
    t_name = ['cp', 'pw']
    
    perf_stats = {}
    for strategy, df in metrics.items():
        inds = ind_tradeoff[strategy]
        best_perf = df.iloc[inds, :]
        best_perf.index = t_name
        describe = df.describe()
        describe = describe.append(best_perf)
        perf_stats[strategy] = describe.copy()
    
    return perf_stats

def get_hv_list(cal_dirs, max_gen):
    
    hv_list = []
    npoints_list = []
    best_list = []
    for i, cal_dir in enumerate(cal_dirs):
        
        hv_gen = []
        npoints = []
        swat_optimal = Optimal(cal_dir)
        
        best_list.append(swat_optimal.get_best())
        
        for ngen in range(1, max_gen+1):
            obj = swat_optimal.get_obj(ngen).to_numpy()
            obj = obj.copy(order='C')
            np_i = obj.shape[0]
            if ngen == 1:
                ref = np.max(obj, axis=0)
            hv = hvwfg.wfg(obj, ref)
            hv_gen.append(hv)
            npoints.append(np_i)
        
        hv_gen = np.array(hv_gen)
        hv_gen = (hv_gen - hv_gen.min())/(hv_gen.max()-hv_gen.min())
        hv_list.append(hv_gen)
        npoints_list.append(npoints)
        
        if i == 0:
            file_list = swat_optimal.obj_ref_list
            for j, file in enumerate(file_list):
                objs = pd.read_csv(cal_dir + '/' + file, header=0)
                constr = objs.copy().to_numpy()[:, 6:]
                temp = constr.copy()
                temp[constr <= 0] = 1
                temp[constr > 0] = 0            
                aux = temp.prod(axis=1)
                
                if sum(aux) > 0:
                    feasible_gen = j + 1
                    break
    
    return hv_list, npoints_list, best_list, feasible_gen

def plot_performance(hv_list, npoints_list, best_list, feasible_gen, output_dir):
    
    fig = plt.figure(figsize=(10,5))
    G = gridspec.GridSpec(2, 2, figure=fig)
    
    # Modified Taylor diagram for Pareto-optimal results
    
    dia = TaylorDiagram(fig=fig, rect=222, label="Reference", srange=(0, 1.8), extend=True)
    dia.add_grid()
    contours = dia.add_contours(levels=5, colors='0.5')
    ls = ['b', 'r', 'g']
    ax0 = plt.subplot(G[1,1])
    ax0.plot([0, 1.8], [1,1], 'k--', label='_nolegend_')
    ax0.plot([1, 1], [0, 1.5], 'k--', label='_nolegend_')
    
    for j, best_s in enumerate(best_list):
        flows = best_s[2]
        flows.columns = list(range(0, len(flows.columns)))
    
        r_list = []
        alpha_list = []
        beta_list = []
        for ind in flows.columns:
            ts = Metric(obs_file, flows[ind])
            ts.get_all()
            
            r_list.append(ts.r)
            alpha_list.append(ts.alpha)
            beta_list.append(ts.beta)
        
        for i, (alpha, corrcoef) in enumerate(zip(alpha_list, r_list)):
            dia.add_sample(alpha, corrcoef, 'k.', markerfacecolor=ls[j], alpha=0.7, markersize=8)
        
        ax0.plot(alpha_list, beta_list, 'k.', markerfacecolor=ls[j], alpha=0.7, markersize=8) 
          
    plt.clabel(contours, inline=1, fontsize=8, fmt='%.2f')
    
    dia.ax.get_shared_x_axes().join(dia.ax, ax0)
    ax0.set_xlim(0,1.8)
    ax0.set_ylim(0,1.5)
    ax0.set_aspect(1.22, adjustable='box', anchor=(0.39,0))
    ax0.set_ylabel(r'$\beta=\mu_s/\mu_o$')
    ax0.set_xlabel(r'$\alpha=\sigma_s/\sigma_o$')
    
    ax0.set_xticks([0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8])
    ax0.set_xticklabels(['0.00','','0.40','','0.80','','1.20','','1.60',''])
    ax0.grid()
    
    ax0.legend(['$\it{Strategy}$ 1', '$\it{Strategy}$ 2', 'Initial population'], bbox_to_anchor=(1,1), loc="upper left")
    ax0.set_title('c)', loc='right')
    
    # Hypervolume and number of solutions plots
    
    ax1 = plt.subplot(G[0,0])
    ax2 = plt.subplot(G[1,0], sharex=ax1)
    plt.setp(ax1.get_xticklabels(), visible=False)
    
    gen = np.arange(1,max_gen+1)
    ls = ['b-','r-']
    ls0 = ['lightblue','lightsalmon']
    for i, hv_i in enumerate(hv_list):
        ax1.plot(gen, hv_i, color=ls0[i], label = '_nolabel_', alpha=0.7)
        txt = '{:d}'.format(i+1)
        ax1.plot(gen, moving_average(hv_i, 10), ls[i], label = '$\it{Strategy}$ '+ txt, alpha=0.7)
        ax2.plot(gen, npoints_list[i], color=ls0[i], label = '$\it{Strategy}$ ' + txt, alpha=0.7)
        ax2.plot(gen, moving_average(npoints_list[i], 10), ls[i], label = '_nolabel_', alpha=0.7)
        if i == 0:
            ax1.plot(feasible_gen, moving_average(hv_i, 10)[feasible_gen-1], 'ko', label='First feasible solution',
                     markerfacecolor='y', markersize=7)
            ax2.plot(feasible_gen, moving_average(npoints_list[i], 10)[feasible_gen-1], 'ko', 
                     label='First feasible solution', markerfacecolor='y', markersize=7)
        
    ax1.set_ylabel('Normalized\n Hypervolume indicator')
    ax2.set_ylabel('Number of near-optimal\n Pareto solutions')
    ax2.set_xlabel('Generation number')
    
    handles, labels = ax1.get_legend_handles_labels()
    order = [0,2,1]
    ax1.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='lower right')
    
    ax1.set_title('a)', loc='right')
    # ax2.set_title('c)', loc='right')
    
    plt.savefig(output_dir + '/MO_Performance.tif', dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    
def get_ERHI_performance(obs_file, best_list, list_ind, groups, output_dir, opt_val=False):
    
    obs = read_obs_file(obs_file)
    obs_hit = Hi(obs)
    eps = np.finfo(float).eps
    
    all_errors = {}
    for group, inds in zip(groups, list_ind):
        for j, best_s in enumerate(best_list[:2]):
            if opt_val:
                flow = best_s
            else:
                flow = best_s[2]
            flow.columns = list(range(0, len(flow.columns)))
            errors = []
            for i in range(0, flow.shape[1]):
                sim = pd.DataFrame(flow.iloc[:, i])
                sim_hit = Hi(sim, obs_file=obs_file)
                index_obs_dict = obs_hit.get_list_indices(inds)
                index_sim_dict = sim_hit.get_list_indices(inds)
                error = [((s+eps) - (o+eps)) / (o+eps) * 100
                          for o, s in zip(list(index_obs_dict.values()), list(index_sim_dict.values()))]
                errors.append(error)
            
            errors = np.array(errors)
            errors = pd.DataFrame(errors, columns=list(index_obs_dict.keys()))
            errors['Strategy'] = np.ones((errors.shape[0], 1))*(j+1)
            if j == 0:
                big_error = errors
            else:
                big_error = pd.concat([big_error, errors], ignore_index=True)
        
        be_melt = big_error.melt(id_vars='Strategy', value_vars=big_error.columns[:-1])    
        
        if opt_val:
            be_melt.to_csv(output_dir + '/' + group + '_val.csv', index=False)
        else:
            be_melt.to_csv(output_dir + '/' + group + '.csv', index=False)
        
        all_errors[group] = big_error
        
    return all_errors

def get_stats(obs_file, best_list, ind_dict, perf_dict, opt_val=False):
    
    obs = read_obs_file(obs_file)
    q_obs = fdc(obs)
    obs_hit = Hi(obs)
    
    perf_strategies = {}
    perf_strategies2 = {}
    
    for j, best_s in enumerate(best_list[:2]):
        
        if opt_val:
            flow = best_s
        else:
            flow = best_s[2]
        
        perf_pareto = []
        perf_pareto2 = []
        column_names = ['FHV', 'FMV', 'FMS', 'FLV', 'RR']
        
        for i in range(0, flow.shape[1]):
            sim = pd.DataFrame(flow.iloc[:, i])
            sim_hit = Hi(sim, obs_file=obs_file)
            sim = sim[obs.index[0]:obs.index[-1]]            
            q_sim = fdc(sim)
            ts = Metric(obs_file, sim)
            
            fdc_pbias = list(fdc_stats(q_obs, q_sim))
            pbias = ts.get_pbias()            
            fdc_pbias.append(pbias)
            
            perf_pareto.append(fdc_pbias)
            
            pm = get_perf(ts, perf_dict)
            hi_perf = get_erhi_mean(obs_hit, sim_hit, ind_dict)
            
            pm.update(hi_perf)
            perf_pareto2.append(pm)
        
        perf_pareto = pd.DataFrame(perf_pareto, columns = column_names)
        perf_pareto2 = pd.DataFrame(perf_pareto2)
        perf_strategies['Strategy_' + str(j + 1)] = perf_pareto.copy()
        perf_strategies2['Strategy_' + str(j + 1)] = perf_pareto2.copy()
    
    return perf_strategies, perf_strategies2

def get_perf(simulation, perf_dict):
    
    perf = {}
    for key, value in perf_dict.items():
        for trans in value:
            name = '{:s} ({:s})'.format(str(key), trans)
            metric = getattr(simulation, 'get_{:s}'.format(str(key).lower()))(trans.lower())
            perf[name] = metric
            
    return perf

def get_erhi_mean(obs_hit, sim_hit, ind_dict):
    
    eps = np.finfo(float).eps
    
    perf = {}
    for key, inds in ind_dict.items():
        index_obs_dict = obs_hit.get_list_indices(inds)
        index_sim_dict = sim_hit.get_list_indices(inds)
        error = [((s+eps) - (o+eps)) / (o+eps) * 100
                  for o, s in zip(list(index_obs_dict.values()), list(index_sim_dict.values()))]
        error = np.abs(np.array(error))
        metric = np.nanmean(error)
        perf[str(key)] = metric
    
    return perf    
        

def plot_tradeoffs(obs_file, obs_file_val, best_list, val_list, metrics, metrics_val, output_dir):
    
    obs = read_obs_file(obs_file)
    
    fig = plt.figure(figsize=(8, 10), constrained_layout=False, tight_layout=True)
    G = gridspec.GridSpec(4, 4, figure=fig, width_ratios=(1,1,0.4,1.6))
    q_obs = fdc(obs)
    col = [['lightblue', 'blue', 'magenta'], ['lightsalmon', 'red']]
    st = ['$\it{Strategy}$ 1', '$\it{Strategy}$ 2']
    tit = ['a)', 'b)', 'c)', 'd)']
    
    ax1 = plt.subplot(G[:2,:3])
    ax1.plot(q_obs, color='black', linestyle='--', label='Observed', linewidth=2, zorder=3)
    
    tradeoff_ind = {}
    
    for j, best_s in enumerate(best_list[:2]):
        # objs = best_s[0].to_numpy()
        objs = best_s[0].to_numpy()
        # params = best_s[1]
        sim = best_s[2]
        # sim = val_list[j]
        sigs = metrics['Strategy_' + str(j + 1)]
        sigs_val = metrics_val['Strategy_' + str(j + 1)]
        
        # best trade-off solution (compromise programming)
        a = compromise(objs.copy(), [0,0,0,0,0,0], p=0, norm=True)[0]
        # best trade-off solution (pseudo-weights)
        weights = np.ones((6,))*1/6
        b, pseudo_weights = get_decision_making("pseudo-weights", weights).do(objs, return_pseudo_weights=True)
        
        tradeoff_ind['Strategy_' + str(j + 1)] = [a, b]
        
        sim = sim[obs.index[0]:obs.index[-1]]
        comp = sim.iloc[:,a]
        pw = sim.iloc[:,b]
        
        q_sim = fdc(sim)
        q_sim_comp = fdc(comp.to_numpy().reshape((-1,1)))
        q_sim_pw = fdc(pw.to_numpy().reshape((-1,1)))
        
        sigs_comp = sigs.iloc[a, :] 
        sigs_pw = sigs.iloc[b, :]
        
        sigs_val_comp = sigs_val.iloc[a, :] 
        sigs_val_pw = sigs_val.iloc[b, :]
        
        ax = plt.subplot(G[j,3])
        ax.plot(q_sim, color='lightblue', linestyle='-', label='Pareto-optimal sols.')
        ax.plot(q_obs, color='black', linestyle='--', label='Observed', linewidth=1.5)
        ax.axvline(x=2, color='gray', linestyle='-', linewidth=1.5, label='_nolegend_')
        ax.axvline(x=20, color='gray', linestyle='-', linewidth=1.5, label='_nolegend_')
        ax.axvline(x=70, color='gray', linestyle='-', linewidth=1.5, label='_nolegend_')
        ax.set_yscale('log')
        # ax.grid(True, which="major", linewidth=1)
        # ax.grid(True, which="minor", linestyle='--')
        ax.set_title(tit[j], loc='right')
        # ax.autoscale(enable=True, axis='x', tight=True)

        
        ax2 = plt.subplot(G[2+j, :2])
        ax2.plot(sigs.transpose(), color='lightblue', linestyle='-', label='Pareto-optimal sols.')
        ax2.plot(sigs_comp, color=col[j][1], linestyle='-.', linewidth=2, 
                 label='Compromise prog. sol. (' + st[j] + ')')
        # ax2.grid(True, which="major", linewidth=1)
        ax2.set_title(tit[j+2], loc='left')
        ax2.set_ylabel('PBIAS (%)')
        ax2.axhline(y=0, color='black', linestyle='--', linewidth=1.5)
        ax2.autoscale(enable=True, axis='x', tight=True)
        
        ax3 = plt.subplot(G[2+j, 2:])
        ax3.plot(sigs_val.transpose(), color='lightblue', linestyle='-', label='Pareto-optimal sols.')
        ax3.plot(sigs_val_comp, color=col[j][1], linestyle='-.', linewidth=2, 
                 label='Compromise prog. sol. (' + st[j] + ')')
        # ax3.grid(True, which="major", linewidth=1)
        ax3.axhline(y=0, color='black', linestyle='--', linewidth=1.5)
        ax3.autoscale(enable=True, axis='x', tight=True)

        # ax2.set_ylim(-120, 120)
        
        if j==0:
            ax.tick_params(labelbottom=False) 
            ax2.tick_params(labelbottom=False)
            ax3.tick_params(labelbottom=False)
            
            ax2.set_ylim(-30,30)
            ax3.set_ylim(-30,30)
            
            ax2.set_title('Calibration')
            ax3.set_title('Validation')
        else:
            ax2.set_ylim(-110,110)
            ax3.set_ylim(-110,110)
        
        ax1.plot(q_sim_comp, color=col[j][1], linestyle='-.', linewidth=2, 
                 label='Compromise prog. sol. (' + st[j] + ')', zorder=2)
        if j == 0:
            ax1.plot(q_sim_pw, color=col[j][2], linestyle='-', linewidth=2, label='Pseudo-weight sol. (' + st[j] + ')',
                     zorder=1)
            ax2.plot(sigs_pw, color=col[j][2], linestyle='-', linewidth=2, 
                 label='Pseudo-weight sol. (' + st[j] + ')')
            ax3.plot(sigs_val_pw, color=col[j][2], linestyle='-', linewidth=2, 
                 label='Pseudo-weight sol. (' + st[j] + ')')
        
    ax1.axvline(x=2, color='gray', linestyle='-', linewidth=1.5, label='_nolegend_')
    ax1.axvline(x=20, color='gray', linestyle='-', linewidth=1.5, label='_nolegend_')
    ax1.axvline(x=70, color='gray', linestyle='-', linewidth=1.5, label='_nolegend_')
    ax1.set_ylabel('Flow ($m^3/s$)')
    ax1.set_xlabel("                                                        Exceedance Pr. (%)")
    # ax1.autoscale(enable=True, axis='x', tight=True)

    
    # ax1.grid(True, which="major", linewidth=1)
    # ax1.grid(True, which="minor", linestyle='--')
    ax1.set_yscale('log')
    ax1.legend(facecolor='white', framealpha=1)
    
    fig.savefig(output_dir + '/FDC.tif', dpi=600, pil_kwargs={"compression": "tiff_lzw"}, bbox_inches='tight')
    
    return tradeoff_ind
    
def plot_PCP(dfs, cols, limits, ynames, n=4, alpha=1): 
    # Based on https://stackoverflow.com/questions/8230638/parallel-coordinates-plot-in-matplotlib
    
    fig, host = plt.subplots(figsize=(10, 5))
    
    df = []
    s = []
    i = 1
    for item in dfs.values():
        temp = item.copy()
        strategy = np.ones((temp.shape[0], 1)) * int(i)
        temp = temp.to_numpy()
        df.append(temp)
        s.append(strategy)
        i = i +1
    
    df = np.vstack(df)
    s = np.vstack(s).astype('int')

    ys = df[:, cols]       
    ymins = np.array([i[0] for i in limits])
    ymaxs = np.array([i[1] for i in limits])

    dys = ymaxs - ymins
    # ymins -= dys * 0.05  # add 5% padding below and above
    # ymaxs += dys * 0.05
    # dys = ymaxs - ymins
    
    # transform all data to be compatible with the main axis
    zs = np.zeros_like(ys)
    zs[:, 0] = ys[:, 0]
    zs[:, 1:] = (ys[:, 1:] - ymins[1:]) / dys[1:] * dys[0] + ymins[0]
    
    axes = [host] + [host.twinx() for i in range(ys.shape[1] - 1)]
    
    for i, ax in enumerate(axes):
        ax.set_ylim(ymins[i], ymaxs[i])
        ax.yaxis.set_major_locator(plt.MaxNLocator(n))
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        if ax != host:
            ax.spines['left'].set_visible(False)
            ax.yaxis.set_ticks_position('right')
            ax.spines["right"].set_position(("axes", i / (ys.shape[1] - 1)))

    host.set_xlim(0, ys.shape[1] - 1)
    host.set_xticks(range(ys.shape[1]))
    host.set_xticklabels(ynames, fontsize=12)
    host.tick_params(axis='x', which='major', pad=7)
    host.spines['right'].set_visible(False)
    host.xaxis.tick_bottom()
    # host.set_title('Parallel Coordinates Plot', fontsize=18)
    
    colors = ['lightblue', 'salmon']
    order = [2, 1]
    for j in range(ys.shape[0]):
        # to just draw straight lines between the axes:
        # host.plot(range(ys.shape[1]), zs[j,:], colors[s[j,0]-1], alpha=alpha, zorder=order[s[j,0]-1])
    
        # create bezier curves
        # for each axis, there will a control vertex at the point itself, one at 1/3rd towards the previous and one
        #   at one third towards the next axis; the first and last axis have one less control vertex
        # x-coordinate of the control vertices: at each integer (for the axes) and two inbetween
        # y-coordinate: repeat every point three times, except the first and last only twice
        verts = list(zip([x for x in np.linspace(0, len(ys) - 1, len(ys) * 3 - 2, endpoint=True)],
                          np.repeat(zs[j, :], 3)[1:-1]))
        # for x,y in verts: host.plot(x, y, 'go') # to show the control points of the beziers
        codes = [Path.MOVETO] + [Path.CURVE4 for _ in range(len(verts) - 1)]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', lw=2, edgecolor=colors[s[j,0]-1], alpha=0.3, 
                                  zorder=order[s[j,0]-1])
        host.add_patch(patch)
        
    plt.tight_layout()

def moving_average(x, w):
    z = np.convolve(x, np.ones(w), 'valid') / w
    z = np.vstack((np.ones((w-1,1)) * z[0], z.reshape((-1,1))))
    return z

def lp_metric(f, z, p):
    
    if len(f) != len(z):
        print('number of dimensions for individual and reference point must be the same')
        return
       
    if p == 0: # Chebyshev metric
        di = [abs(fi-zi) for fi, zi in zip(f, z)]
        d = max(di)
    else: # lp metric
        di = [abs(fi-zi)**p for fi, zi in zip(f, z)]
        d = sum(di)**(1/p)
    
    return d

def compromise(obj, z, p, norm=True):
    
    if type(obj) == 'list':
        obj = np.array(obj)
    
    if norm:
        obj += -(np.min(obj, axis=0))
        obj /= np.max(obj, axis=0)
    
    d = [lp_metric(x.tolist(),z,p) for x in obj]
    best = [i for i, x in enumerate(d) if x == np.min(d)]
    
    return best

def fdc_interp(q, xp):
    q = q.sort_index()
    
    x = np.array(q.index)
    y = np.array(q.values)[:,0]
    yp = np.interp(xp, x, y)
    
    return yp

def bias_fhv(obs, sim, interval):
    pr = np.arange(interval[0], interval[1], 0.01)    
    obs = fdc_interp(obs, pr)
    sim = fdc_interp(sim, pr)
    
    res = sim - obs
    pbias = res.sum() / obs.sum() * 100   
    
    return pbias

def bias_flv(obs, sim, interval):
    min_obs = obs.to_numpy().min()
    min_sim = sim.to_numpy().min()
    min_q = np.min((min_obs, min_sim))
    
    pr = np.arange(interval[0], interval[1], 0.01)    
    obs = fdc_interp(obs, pr)
    sim = fdc_interp(sim, pr)
    
    vsim = (np.log(obs) - np.log(min_q)).sum()
    vobs = (np.log(sim) - np.log(min_q)).sum()
    pbias = -1 * (vsim - vobs) / vobs * 100
    
    return pbias

def bias_ms(obs, sim, interval):
    obs = fdc_interp(obs, interval)
    sim = fdc_interp(sim, interval)
    
    s_sim = np.log(sim[0]) - np.log(sim[1])
    s_obs = np.log(obs[0]) - np.log(obs[1])
    
    pbias = (s_sim - s_obs) / s_obs * 100
    
    return pbias

def fdc_stats(obs, sim):
    peak = [0, 2]
    high = [2, 20]
    mid = [20, 70]
    low = [70, 100]
    
    pbias_peak = bias_fhv(obs, sim, peak)
    pbias_high = bias_fhv(obs, sim, high)
    pbias_midSlope = bias_ms(obs, sim, mid)
    pbias_low = bias_flv(obs, sim, low)
    
    return pbias_peak, pbias_high, pbias_midSlope, pbias_low


class TaylorDiagram(object): # Taylor diagram code adapted from https://gist.github.com/ycopin/3342888

    """
    Taylor diagram.
    Plot model standard deviation and correlation to reference (data)
    sample in a single-quadrant polar plot, with r=alpha (ratio of simulated and observed stdev) and
    theta=arccos(correlation).
    """

    def __init__(self, fig=None, rect=111, refstd=1, label='_nolegend_', srange=(0, 1.5), extend=False):
        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.
        Parameters:
        * refstd: reference alpha to be compared to
        * fig: input Figure or None
        * rect: subplot definition
        * label: reference label
        * srange: stddev axis extension, in units of *refstd*
        * extend: extend diagram to negative correlations
        """

        from matplotlib.projections import PolarAxes
        import mpl_toolkits.axisartist.floating_axes as FA
        import mpl_toolkits.axisartist.grid_finder as GF

        self.refstd = refstd            # Reference standard deviation

        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.array([0, 0.2, 0.4, 0.6, 0.8, 0.9, 1])
        self.tmin = 0
        if extend:
            # Diagram extended to negative correlations
            self.tmax = np.arccos(-0.4)
            rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
        else:
            # Diagram limited to positive correlations
            self.tmax = np.pi/2
        tlocs = np.arccos(rlocs)        # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)    # Positions
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0] * self.refstd
        self.smax = srange[1] * self.refstd

        ghelper = FA.GridHelperCurveLinear(
            tr,
            extremes=(self.tmin, self.tmax, self.smin, self.smax),
            grid_locator1=gl1, tick_formatter1=tf1)
        
        if fig is None:
            fig = plt.figure()

        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")   # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Correlation $\it{r}$")

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].toggle(ticklabels=False, label=False)
        ax.axis["left"].label.set_text(r'$\alpha=\sigma_s/\sigma_o$')

        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=False, label=False)
        ax.axis["right"].label.set_text(r'$\alpha=\sigma_s/\sigma_o$')
        ax.axis["right"].major_ticklabels.set_axis_direction(
            "bottom" if extend else "left")
        
        ax.set_anchor('SW')
        ax.set_title('b)', loc='right')

        if self.smin:
            ax.axis["bottom"].toggle(ticklabels=False, label=False)
        else:
            ax.axis["bottom"].set_visible(False)          # Unused

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        # Add reference point and stddev contour
        l, = self.ax.plot([0], self.refstd, 'ko',
                          ls='', ms=3, label=label)
        t = np.linspace(self.tmin, self.tmax)
        r = np.zeros_like(t) + self.refstd
        self.ax.plot(t, r, 'k--', label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]

    def add_sample(self, alp, corrcoef, *args, **kwargs):
        """
        Add sample (*alpha*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        `Figure.plot` command.
        """

        l, = self.ax.plot(np.arccos(corrcoef), alp,
                          *args, **kwargs)  # (theta, radius)
        self.samplePoints.append(l)

        return l

    def add_grid(self, *args, **kwargs):
        """Add a grid."""

        self._ax.grid(*args, **kwargs)
        
    def add_contours(self, levels=5, **kwargs):
        """
        Add constant centered RMS difference contours, defined by *levels*.
        """

        rs, ts = np.meshgrid(np.linspace(self.smin, self.smax),
                             np.linspace(self.tmin, self.tmax))
        # Compute centered RMS difference
        rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)

        return contours
  
    
if __name__ ==  '__main__':
    
    obs_file = os.path.abspath('../resources/Observed/Honeyoy.csv')
    obs_file_val = os.path.abspath('../resources/Observed/Honeyoy_8090.csv')

    list_ind = [['ma12', 'ma13', 'ma14', 'ma15', 'ma16', 'ma17', 'ma18', 'ma19', 'ma20', 'ma21', 'ma22', 'ma23'],
                ['dl1', 'dl2', 'dl3', 'dl4', 'dl5', 'dh1', 'dh2', 'dh3', 'dh4', 'dh5', 'ml17'], ['tl1', 'th1',
                'fl1', 'dl16', 'fh1', 'dh15', 'ra1', 'ra3', 'ra8'], 
                ['mag1', 'mag2', 'mag3', 'mag4', 'mag5', 'mag6', 'mag7']]
    
    groups = ['IHA1', 'IHA2', 'IHA3_4_5', 'MAG7']
    
    list_ind2 = [['ma24', 'ma25', 'ma26', 'ma27', 'ma28', 'ma29', 'ma30', 'ma31', 'ma32', 'ma33', 'ma34', 'ma35'],
                ['dl6', 'dl7', 'dl8', 'dl9', 'dl10', 'dh6', 'dh7', 'dh8', 'dh9', 'dh10', 'ml18'],
                ['tl2', 'th2','fl2', 'dl17', 'fh2', 'dh16', 'ra2', 'ra4', 'ra9']]
    
    groups2 = ['IHA1', 'IHA2', 'IHA3_4_5']
    
    ind_dict = {'IHA1': ['ma12', 'ma13', 'ma14', 'ma15', 'ma16', 'ma17',
                         'ma18', 'ma19', 'ma20', 'ma21', 'ma22', 'ma23'],
                'IHA2': ['dl1', 'dl2', 'dl3', 'dl4', 'dl5', 'dh1',
                         'dh2', 'dh3', 'dh4', 'dh5', 'ml17'],
                'IHA3': ['tl1', 'th1'],
                'IHA4': ['fl1', 'dl16', 'fh1', 'dh15'],
                'IHA5': ['ra1', 'ra3', 'ra8'],
                'MAG7': ['mag1', 'mag2', 'mag3', 'mag4', 'mag5', 'mag6', 'mag7']}
    
    ind_dict2 = {'IHA1': ['ma24', 'ma25', 'ma26', 'ma27', 'ma28', 'ma29',
                          'ma30', 'ma31', 'ma32', 'ma33', 'ma34', 'ma35'],
                'IHA2': ['dl6', 'dl7', 'dl8', 'dl9', 'dl10', 'dh6',
                         'dh7', 'dh8', 'dh9', 'dh10', 'ml18'],
                'IHA3': ['tl2', 'th2'],
                'IHA4': ['fl2', 'dl17', 'fh2', 'dh16'],
                'IHA5': ['ra2', 'ra4', 'ra9']}
    
    perf_dict = {'KGE': ['none', 'inverse'],
                 'R2': ['none', 'sqrt', 'inverse'],
                 'IoA': ['none', 'rel'],
                 'NSE': ['none']}
    
    val_file = os.path.abspath('../output/validation.pickle')
    
    # Processing calibration results
    cal_dirs = ['../output/Problem_0', '../output/Problem_3']
    max_gen = 1000
    
    opt = [True, True, True]
    
    output_cal, output_val, metrics_cal,  metrics_val, errors_cal, errors_val = main(obs_file, obs_file_val, 
                                                                                     list_ind, groups, ind_dict, 
                                                                                     perf_dict, cal_dirs, max_gen, 
                                                                                     val_file, opt)