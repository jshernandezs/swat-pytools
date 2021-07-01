clc, clear, close all

%% Settings
rng('default')

nprocessors = 4;
Extra.n_model_params = 15;
Extra.path_to_python = '~/anaconda3/envs/swatpy/bin/python';
Extra.lambda = 0.2; % Box-cox transformation parameter
Extra.psi = 0.8; % Autocorrelation error model parameter
Extra.wrapper = 'swat_wrapper.py';

%% Define input x variable

nsim = 10;

x_min = [0, -0.25, 0, 0, 0, 0, 0, 0, 0.02, 0, 0, 0, 0, -0.25, 1, 0.5];
x_max = [1, 0.25, 100, 1, 1, 500, 1, 5000, 0.2, 1000, 1, 0.3, 500, 0.25, 24, 10];

x = repmat(x_min, nsim, 1) + lhsdesign(nsim, length(x_min)).*repmat(x_max - x_min, nsim, 1);

%% Run Monte Carlo Simulation

parpool(nprocessors);

logL = zeros(nsim, 1);
parfor i=1:nsim
    logL(i) = loglikelihood(x(i, :), Extra);
end

delete(gcp('nocreate'))