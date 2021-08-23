clc, clear, close all

%% Settings

Extra.n_model_params = 15;
Extra.path_to_python = '~/anaconda3/envs/swatpy/bin/python';
Extra.lambda = 0.5; % Box-cox transformation parameter
Extra.psi = 0.5; % Autocorrelation error model parameter

%% Define input x variable
x = [0.1, -0.21, 1.67, 0.70, 0.0059, 6.11, 0.83, 438, 0.16, 438,...
     0.50, 0.12, 6.45, -0.21, 1.10, 1];

%% Compute log likelihood (McInerney et al., 2018) 
logL = loglikelihood(x, Extra);