%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%                                                                                                          %
%% DDDDDD   RRRRRR   EEEEEEE   AAAA   MMM   MMM PPPPPPP    AAAA   CCCCCCC KKK  KKK   AAAA   GGGGGGG EEEEEEE %
%% DDDDDDD  RRRRRRR  EEEEEEE  AAAAAA  MMM   MMM PPPPPPPP  AAAAAA  CCCCCC  KKK  KKK  AAAAAA  GGG     EEEEEEE % 
%% DDD  DDD RRR RRR  EEE     AAA  AAA MMMM MMMM PPP  PPP AAA  AAA CCC     KKK KKK  AAA  AAA GGG     EEE     % 
%% DDD  DDD RRR RRR  EEEEEE  AAA  AAA MMMMMMMMM PPP  PPP AAA  AAA CCC     KKKKKK   AAA  AAA GGG GGG EEEEEE  %
%% DDD  DDD RRRRRRR  EEEEEE  AAAAAAAA MMMM MMMM PPPPPPPP AAAAAAAA CCC     KKKKKK   AAAAAAAA GGG GGG EEEEEE  %
%% DDD  DDD RRR RRR  EEE     AAAAAAAA MMM   MMM PPPPPPP  AAAAAAAA CCC     KKK KKK  AAAAAAAA GGG GGG EEE     %
%% DDDDDDD  RRR  RRR EEEEEEE AAA  AAA MMM   MMM PPP      AAA  AAA CCCCCC  KKK  KKK AAA  AAA GGGGGGG EEEEEEE %
%% DDDDDD   RRR  RRR EEEEEEE AAA  AAA MMM   MMM PPP      AAA  AAA CCCCCCC KKK  KKK AAA  AAA GGGGGGG EEEEEEE %
%%                                                                                                          %
%% DREAM PACKAGE: DiffeRential Evolution Adaptive Metropolis algorithm with continuous and/or               %
%% discrete variables and sampling from an archive of current or past states. This package implements in    %
%% in one code the DREAM, DREAM_D, DREAM_ZS and DREAM_ZS algorithms. These four members of the DREAM family %
%% used to have their own separate toolbox but are now merged (release January 2018) into one package.      %
%%                                                                                                          %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%                                                                                                          %
%% SYNOPSIS:                                                                                                %
%%  [chain,output,FX,logL] = DREAM_Package(method,Func_name,DREAMPar,Par_info)                              %
%%  [chain,output,FX,logL] = DREAM_Package(method,Func_name,DREAMPar,Par_info,Meas_info)                    %
%%  [chain,output,FX,logL] = DREAM_Package(method,Func_name,DREAMPar,Par_info,Meas_info,options)            %
%%  [chain,output,FX,logL] = DREAM_Package(method,Func_name,DREAMPar,Par_info,Meas_info,options,plug)       %
%%                                                                                                          %
%% INPUT ARGUMENTS:                                                                                         %
%%  method      REQUIRED: Name (string) of MCMC algorithm ('DREAM'/'DREAM_ZS'/'DREAM_D'/'DREAM_DZS')        %
%%  Func_name   REQUIRED: Function name (string) that returns simulation/(log)-likelihood/(log)-density     %
%%  DREAMPar    REQUIRED: Structure with algorithmic settings of MCMC algorithm                             %
%%   .d             Dimensionality (# variables) target distribution                                        %
%%   .N             Number of Markov chains                                                                 %
%%   .T             Number of generations                                                                   %
%%   .lik           Choice of likelihood function                                                           %
%%   .nCR           Number of crossover values                      (default: 3)                            %
%%   .delta         Number chain pairs for proposal                 (default: 3)                            %
%%   .lambda        Random error for ergodicity                     (default: 0.05)                         %
%%   .zeta          Randomization                                   (default: 1e-12)                         %
%%   .p_unit_gamma  Selection probability unit jumprate (gamma)     (default: 0.2)                          %
%%   .adapt_pCR     Adapt selection prob. crossover?                (default: 'yes')                        %
%%   .thinning      Each Tth sample is stored                       (default: 1)                            %
%%   .GLUE          GLUE likelihood parameter                       (default: 10)                           %
%%   .beta0         Scaling factor of built-in jump rate            (default: 1)                            %
%%   .outlier       Outlier chain detection test                    (default: 'iqr') -> DREAM/DREAM_D       %
%%   .psnooker      Selection probability of snooker jump           (default: 0.1)   -> DREAM_ZS/DREAM_DZS  %
%%   .m0            Initial size of external archive, Z             (default: 10*d)  -> DREAM_ZS/DREAM_DZS  %
%%   .k             Growth rate of external archive                 (default: 10)    -> DREAM_ZS/DREAM_DZS  %
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
%%   .steps         d-vector with # intervals of each parameter                      -> DREAM_D/DREAM_DZS   %
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
%% OUTPUT ARGUMENTS:                                                                                        %
%%  chain       Three-dimensional array with N chain trajectories, log-prior and log-likelihood values      %
%%  output      Structure that summarizes algorithmic performance                                           %
%%   .R_stat    Matrix with evolution univariate \hat{R} convergence diagnostic                             %
%%   .MR_stat   Matrix with evolution multivariate \hat{R} convergence diagnostic                           %
%%   .AR        Matrix with evolution of acceptance rate (%)                                                %
%%   .CR        Matrix with evoluton of crossover selection probabilities                                   %
%%   .outlier   Matrix with generation number and index of outlier chains                                   %
%%   .RunTime   Total required CPU time in seconds                                                          %
%%  FX          Matrix with model simulations (listed as rows)                                              %
%%  logL        Matrix with log-likelihood values of sampled chains (column wise)                           %
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
%%                                                                                                          %
%% FURTHER CHECKING:                                                                                        %
%%  Website:  http://faculty.sites.uci.edu/jasper                                                           %
%%  Papers: http://faculty.sites.uci.edu/jasper/publications/                                               %
%%  Google Scholar: https://scholar.google.com/citations?user=zkNXecUAAAAJ&hl=nl                            %                                                                                                        %
%% 
%% HOW TO INSTALL TOOLBOX AND EXECUTE EXAMPLES:                                                             %
%%  Step 1: Go to main DREAM PACKAGE directory                                                              %
%%  Step 2: Type "addpath(pwd,[pwd '/postprocessing'],[pwd '/diagnostics'],[pwd '/gamesampling']);"         %
%%  Step 3: Go to an example directory; say example_1: Type "cd example_1"                                  %
%%  Step 4: Execute this example by typing in command prompt: "example_1"                                   %
%%                                                                                                          %
%%  After DREAM PACKAGE has terminated you get visual output (if options.print = 'yes');                    %
%%  You can create a 2d matrix from 3D-chain array by typing "ParSet = genparset(chain);                    %
%%  This gives you all sampled solutions (apply burn-in before you create histograms of "Pars")             %
%%                                                                                                          %
%% HOW TO USE GAME SAMPLING?                                                                                %
%%  Step 1: Get posterior samples, type, say "Xp = Pars(end-15000:end,1:DREAMPar.d+2);"                     % 
%%  Step 2: Now calculate marginal likelihood using "Z = GAME_sampling(Xp,method,DREAMPar,Func_name);"      %
%%          where method is any option from {'ris'/'is'/'ob'/'gb'}                                          %
%%          'ris' : reciprocal importance sampling                                                          %
%%          'is'  : importance sampling                                                                     %
%%          'ob'  : optimal bridge sampling                                                                 %
%%          'gb'  : geometric bridge sampling                                                               % 
%%  Type in MATLAB prompt "edit readme.txt" for information about how to use GAME_sampling                  %                                                               %
%%                                                                                                          %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clc, clear, close all
%% Go to main DREAM PACKAGE directory 
addpath(pwd,[pwd '/postprocessing'],[pwd '/diagnostics'],[pwd '/gamesampling']);
%% Now go to problem directory;
cd swat_dream
%% Now execute the problem: 
run_DREAM;
%% After DREAM PACKAGE has terminated you can create a 2d matrix from 3D-chain array as follows
ParSet = GenParSet(chain);
%% And you can initiate GAME_sampling via help GAME_sampling