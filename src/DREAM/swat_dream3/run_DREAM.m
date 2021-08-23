rng('default')

%% Problem settings defined by user
DREAMPar.d = 16;                           % Dimension of the problem
DREAMPar.lik = 2;                          % Custom log-likelihood function (McInerney et al., 2018) 

%% Provide information parameter space and initial sampling
Par_info.initial = 'logprior';                                                  % Prior distribution - multiobjective opt. results
Par_info.logprior = @(x,a,b) logmvkpdf(x,a,b);
Par_info.a = load('po_prior.csv');
Par_info.b = [0.0065,0.0066,1.89,0.017,0.018,0.27,0.021,0.83,0.0041,1.23,0.017,0.00097,0.059,0.030,9.7e-05,0.0043];

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
        DREAMPar.T = 200;                                                     % Number of generations        
end

%% Optional settings
options.modout = 'yes';                                                         % Return model (function) simulations of samples (yes/no)?
options.parallel = 'yes';                                                       % Run each chain on a different core
options.save = 'yes';                                                           % Save workspace DREAM during run

%% Plugin settings

plugin.n_model_params = 15;
plugin.path_to_python = '~/anaconda3/envs/swatpy/bin/python';
plugin.lambda = 0.2;                                                            % Box-cox transformation parameter
plugin.psi = 0.8;                                                               % Autocorrelation error model parameter
plugin.wrapper = 'swat_wrapper.py';                                             % python file containing model to run

%% Run DREAM package
if strcmp(method, 'mtdream_zs')
    [chain,output,FX,Z] = MTDREAM_ZS(Func_name,DREAMPar,Par_info,[],options,plugin);
else
    [chain,output,FX,Z,logL] = DREAM_package(method,Func_name,DREAMPar,Par_info,[],options,plugin);
end