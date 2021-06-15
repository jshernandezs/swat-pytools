%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%                                                                                                          %
%%          MMM   MMM TTTTTTTTT       DDDDDD   RRRRRR   EEEEEEE   AAAA   MMM   MMM                          %
%%          MMM   MMM TTTTTTTTT       DDDDDDD  RRRRRRR  EEEEEEE  AAAAAA  MMM   MMM                          % 
%%          MMMM MMMM TT TTT TT       DDD  DDD RRR RRR  EEE     AAA  AAA MMMM MMMM                          % 
%%          MMMMMMMMM    TTT     ---  DDD  DDD RRR RRR  EEEEEE  AAA  AAA MMMMMMMMM ZZZZZZ SSSSSS            %
%%          MMMM MMMM    TTT     ---  DDD  DDD RRRRRRR  EEEEEE  AAAAAAAA MMMM MMMM    ZZZ SSS               %
%%          MMM   MMM    TTT          DDD  DDD RRR RRR  EEE     AAAAAAAA MMM   MMM   ZZZ  SSSSSS            %
%%          MMM   MMM    TTT          DDDDDDD  RRR  RRR EEEEEEE AAA  AAA MMM   MMM ZZZ       SSS            %
%%          MMM   MMM    TTT          DDDDDD   RRR  RRR EEEEEEE AAA  AAA MMM   MMM ZZZZZZ SSSSSS            %
%%                                                                                                          %
%% MT-DREAM(ZS) is a Markov Chain Monte Carlo algorithm that uses multi-try sampling from an archive        %
%% of past states to speed-up convergence of multiple chains to the target distribution. The algorithm is   % 
%% an extension of the DREAM(ZS) algorithm to multi-try sampling. Diversity of candidate points is enhanced %
%% by the use of an snooker proposal distribution.                                                          %        
%%                                                                                                          %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%                                                                                                          %
%% SYNOPSIS:                                                                                                %
%%  [chain,output,FX,Z] = MTDREAM_ZS(Func_name,DREAMPar,Par_info)                                           %
%%  [chain,output,FX,Z] = MTDREAM_ZS(Func_name,DREAMPar,Par_info,Meas_info)                                 %
%%  [chain,output,FX,Z] = MTDREAM_ZS(Func_name,DREAMPar,Par_info,Meas_info,options)                         %
%%  [chain,output,FX,Z] = MTDREAM_ZS(Func_name,DREAMPar,Par_info,Meas_info,options,plug)                    %
%%                                                                                                          %
%% INPUT ARGUMENTS:                                                                                         %
%%  Func_name   REQUIRED: Function name (string) that returns simulation/(log)-likelihood/(log)-density     %
%%  DREAMPar    REQUIRED: Structure with algorithmic settings of MCMC algorithm                             %
%%   .d             Dimensionality (# variables) target distribution                                        %
%%   .N             Number of Markov chains                         (default: 3)                            %
%%   .T             Number of generations                                                                   %
%%   .lik           Choice of likelihood function                                                           %
%%   .nCR           Number of crossover values                      (default: 3)                            %
%%   .delta         Number of chain pairs for proposal creation     (default: 3)                            %
%%   .lambda        Random error for ergodicity                     (default: 0.05)                         %
%%   .zeta          Randomization                                   (default: 0.05)                         %
%%   .p_unit_gamma  Selection probability unit jumprate (gamma)     (default: 0.2)                          %
%%   .adapt_pCR     Adapt selection prob. crossover?                (default: 'yes')                        %
%%   .thinning      Each Tth sample is stored                       (default: 1)                            %
%%   .GLUE          GLUE likelihood parameter                       (default: 10)                           %
%%   .beta0         Scaling factor of built-in jump rate            (default: 1)                            %
%%   .psnooker      Selection probability of snooker jump           (default: 0.1)                          %
%%   .m0            Initial size of external archive, Z             (default: 10*d)                         %
%%   .k             Growth rate of external archive                 (default: 10)                           %
%%   .mt            Number of multi-try proposals                   (default: 5)                            %
%%  Par_info    REQUIRED: Structure with parameter ranges, initial/prior distribution and bound handling    %
%%   .initial       Method to draw initial chain states                                                     %
%%   .min           d-vector with minimum values of the parameters  (default: -inf*ones(1,DREAMPar.d) )     %
%%   .max           d-vector with maximum values of the parameters  (default: inf*ones(1,DREAMPar.d) )      %   
%%   .mu            Mean of initial sampling distributon            (.initial = 'normal')                   %
%%   .cov           Covariance of initial sampling distribution     (.initial = 'normal')                   %
%%   .x0            N x d matrix with initial chain states          (.initial = 'user')                     %
%%   .prior         Prior distribution (see manual)                 (.initial = 'prior')                    %
%%   .logprior      Prior distribution (see erratum document)       (.initial = 'logprior')                 %
%%   .boundhandling Boundary method ('reflect'/'bound'/'fold')      (default: 'none')                       % 
%%  Meas_info   OPTIONAL: Structure with measurement information (for inference)                            %
%%   .Y             Nx1 vector against which model output is compared                                       %
%%   .Sigma         Nx1 vector with corresponding measurement error                                         %
%%   .S             Scalar/vector with summary metrics                                                      %
%%  options     OPTIONAL: Structure with additional settings/options                                        %
%%   .rho           ABC distance function ( inline format )         (default: inline('X-Y'))                %
%%   .epsilon       ABC epsilon value ( scalar or vector )          (default: 0.025)                        %
%%   .DB            Diagnostic Bayes?                               (default: 'no')                         %
%%   .parallel      Multi-core computation chains?                  (default: 'yes')                        %
%%   .IO            If parallel, IO writing model?                  (default: 'no')                         %
%%   .modout        Return model (function) simulations?            (default: 'no')                         %
%%   .save          Save DREAM output during the run?               (default: 'no')                         %
%%   .restart       Restart run? (only with "save")                 (default: 'no')                         %
%%   .diagnostics   Compute within-chain diagnostics?               (default: 'yes')                        %
%%   .print         Output writing to screen (tables/figures)       (default: 'yes')                        %
%%   .burnin        Burn-in % chain for convergence diagnostics     (default: 50)                           %
%%  plug        OPTIONAL: Second input argument Func_name. Class determined by user                         %
%%                                                                                                          %
%%                                                                                                          %
%% THE DIFFERENT ALGORITHMS HAVE BEEN DESCRIBED IN:                                                         %
%%   Vrugt, J.A., Markov chain Monte Carlo simulation using the DREAM software package: Theory, concepts,   %
%%       and MATLAB Implementation, Environmental Modelling & Software, 75, 273-316,                        %
%%       doi:10.1016/j.envsoft.2015.08.013.                                                                 %
%%   Sadegh, M., and J.A. Vrugt, Approximate Bayesian computation using Markov chain Monte Carlo            %
%%       simulation: DREAM_(ABC), Water Resources Research, doi:10.1002/2014WR015386, 2014                  %
%%   Vrugt, J.A., and M. Sadegh, Toward diagnostic model calibration and evaluation: Approximate Bayesian   %
%%       computation, Water Resources Research, 49, 4335�4345, doi:10.1002/wrcr.20354, 2013.                %
%%   Laloy, E., and J.A. Vrugt, High-dimensional posterior exploration of hydrologic models using multiple- %
%%       try DREAM_(ZS) and high-performance computing, Water Resources Research, 48, W01526,               %
%%       doi:10.1029/2011WR010608, 2012.                                                                    %
%%   Vrugt, J.A., and C.J.F. ter Braak, DREAM_(D): An adaptive Markov chain Monte Carlo simulation          %
%%      algorithm to solve discrete, noncontinuous, and combinatorial posterior parameter estimation        %
%%      problems, Hydrology and Earth System Sciences, 15, 3701-3713, doi:10.5194/hess-15-3701-2011, 2011.  %
%%   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson, Equifinality of formal (DREAM) and       %
%%       informal (GLUE) Bayesian approaches in hydrologic modeling?, Stochastic Environmental Research and %
%%       Risk Assessment, 23(7), 1011-1026, doi:10.1007/s00477-008-0274-y, 2009.                            %
%%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman, Accelerating     %
%%       Markov chain Monte Carlo simulation by differential evolution with self-adaptive randomized        %
%%       subspace sampling, International Journal of Nonlinear Sciences and Numerical Simulation, 10(3),    %
%%       271-288, 2009.                                                                                     %
%%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson, Treatment of input           %
%%      uncertainty in hydrologic modeling: Doing hydrology backward with Markov chain Monte Carlo          %
%%      simulation, Water Resources Research, 44, W00B09, doi:10.1029/2007WR006720, 2008.                   %
%%   Ter Braak, C.J.F., and J.A. Vrugt, Differential Evolution Markov Chain with snooker updater and fewer  %
%%      chains, Statistics and Computing, 10.1007/s11222-008-9104-9, 2008.                                  %
%%   Ter Braak, C.J.F., A Markov Chain Monte Carlo version of the genetic algorithm Differential Evolution: %
%%      easy Bayesian computing for real parameter spaces, Statistics and Computing, 16, 239 - 249,         %
%%      doi:10.1007/s11222-006-8769-1, 2006.                                                                %
%%                                                                                                          %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

%% MT-DREAM_(ZS) developed by Jasper A. Vrugt and Eric Laloy                                     
%% Check:  http://faculty.sites.uci.edu/jasper
%% Papers: http://faculty.sites.uci.edu/jasper/publications/
%% Google Scholar: https://scholar.google.com/citations?user=zkNXecUAAAAJ&hl=nl

%% Different test examples
%% example 1: n-dimensional banana shaped Gaussian distribution
%% example 2: n-dimensional Gaussian distribution
%% example 3: n-dimensional multimodal mixture distribution
%% example 4: real-world example rainfall-runoff (hymod coded in MATLAB)
%% example 5: real-world example rainfall-runoff (hymod coded in FORTRAN)
%% example 6: rainfall-runoff with generalized log-likelihood function
%% example 7: HYDRUS-1D soil hydraulic model: multiplicative prior
%% example 8: 1D mixture distribution: Approximate Bayesian Computation
%% example 9: Rainfall-runoff model with spectral likelihood function
%% example 10: multimodel mixture distrbibution: multivariate prior
%% example 11: multivariate student t distribution
%% example 12: pedometrics problem involving variogram fitting
%% example 13: Nash-Cascade hydrograph
%% example 15: ABC inference using 10 bivariate normal distributions
%% example 16: Hydrogeophysical inference
%% example 17: rainfall - runoff with inference discharge data error
%% example 18: informal likelihood (GLUE): Lotka-Volterra model

%% ------------------------------------------------------------------------
%% Most of the examples herein are used to illustrate how the code works.
%% Yet, for some of the problems it is much better to run DREAM package as
%% this will provide many more posterior samples and hence better posterior
%% estimates. The multi-try sampling per definition stores far fewer samples
%% for a given number of function evaluations, and hence this affects the
%% posterior estimates. Hence please use DREAM package for simple problems
%% For distributed problems storage is sufficient as the the number of 
%% function evaluations will typically be very large.
%% ------------------------------------------------------------------------
clc, clear, close all
%% Go to main MT-DREAM_ZS directory 
addpath(pwd,[pwd '/postprocessing'],[pwd '/diagnostics'],[pwd '/gamesampling']);
%% Now go to example directory; say example_1
cd ../
cd DREAM/swat_dream
%% Now execute this example by typing in command prompt: "example_1" 
run_DREAM;
%% After MT-DREAM_ZS is done you can create a 2d matrix from 3D-chain array
%% ParSet = GenParSet(chain);
%% And you can initiate GAME_sampling via help GAME_sampling