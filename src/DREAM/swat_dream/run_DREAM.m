rng('default')

%% Problem settings defined by user
DREAMPar.d = 16;                           % Dimension of the problem
DREAMPar.lik = 2;                          % Custom log-likelihood function (McInerney et al., 2018) 

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';                                                     % Latin hypercube sampling
Par_info.boundhandling = 'reflect';                                             % Explicit boundary handling
Par_info.min = [0,-0.25,0,0,0,0,0,0,0.02,0,0,0,0,-0.25,1,0.5];                  % If 'latin', min values
Par_info.max = [1,0.25,100,1,1,500,1,5000,0.2,1000,1,0.3,500,0.25,24,5];        % If 'latin', max values

%% Define name of function (.m file) for posterior exploration
Func_name = 'loglikelihood';

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs'}
method = 'mtdream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                                                        % Number of Markov chains
        DREAMPar.T = 20000;                                                     % Number of generations
    case {'dream_zs','dream_dzs'}
        DREAMPar.N = 3;                                                         % Number of Markov chains
        DREAMPar.T = 20000;                                                     % Number of generations
    case {'mtdream_zs'}
        DREAMPar.N = 3;                                                         % Number of Markov chains
        DREAMPar.T = 10000;                                                     % Number of generations        
end

%% Optional settings
options.modout = 'yes';                                                         % Return model (function) simulations of samples (yes/no)?
options.parallel = 'yes';                                                       % Run each chain on a different core
options.save = 'yes';                                                           % Save workspace DREAM during run

%% Plugin settings

plugin.n_model_params = 15;
plugin.path_to_python = '/work/herna505/anaconda3/envs/swatpy/bin/python';
plugin.lambda = 0.2;                                                            % Box-cox transformation parameter
plugin.psi = 0.8;                                                               % Autocorrelation error model parameter

%% Run DREAM package
if strcmp(method, 'mtdream_zs')
    [chain,output,fx,Z] = MTDREAM_ZS(Func_name,DREAMPar,Par_info,[],options,plugin);
else
    [chain,output,FX,Z,logL] = DREAM_package(method,Func_name,DREAMPar,Par_info,[],options,plugin);
end