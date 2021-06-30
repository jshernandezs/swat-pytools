% Script for sampling quantities of interest
clc, clear, close all
rng('default')

cd ..
addpath(pwd);
cd swat_dream

%% settings
filename = 'MTDREAM_Results/cal_results_noninf.mat';
burnin = 80;
nprocessors = 20;

%% get parameter sets
load(filename)

d = DREAMPar.d;
n_model_param = plugin.n_model_params;
% plugin.path_to_python = '~/anaconda3/envs/swatpy/bin/python';

ParSet = GenParSet(chain);
Nset = length(ParSet);

ParEval = ParSet(floor(burnin * Nset / 100) + 1:end, 1:d);
Neval = length(ParEval);

%% run predictions (parameter and total error)

% run the first model (we need it to get number of time steps)
[model, pred, obs] = pred_sim(ParEval(1, :), plugin);
Nt = length(obs);

parpool(nprocessors);

models = zeros(Neval, Nt);
preds = models;

models(1,:) = model;
preds(1,:) = pred;

parfor i=2:Neval
    [models(i,:), preds(i,:)] = pred_sim(ParEval(i, :), plugin);
end

delete(gcp('nocreate'))

save('predictions','preds','models','obs')