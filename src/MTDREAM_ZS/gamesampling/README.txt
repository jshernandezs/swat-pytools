%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%                                                                                                          %
%% This function calculates the marginal likelihood from a collection of target samples                     %
%%                                                                                                          %
%% SYNOPSIS:                                                                                                %
%%   Z = GAME_sampling(Xp,method,DREAMPar,Func_name);                                                       %
%%   Z = GAME_sampling(Xp,method,DREAMPar,Func_name,GAME_options);                                          %
%%   Z = GAME_sampling(Xp,method,DREAMPar,Func_name,GAME_options,Par_info);                                 %
%%   Z = GAME_sampling(Xp,method,DREAMPar,Func_name,GAME_options,Par_info,Meas_info);                       %
%%   Z = GAME_sampling(Xp,method,DREAMPar,Func_name,GAME_options,Par_info,Meas_info,options);               %
%%   Z = GAME_sampling(Xp,method,DREAMPar,Func_name,GAME_options,Par_info,Meas_info,options,plug);          %
%%                                                                                                          %
%% INPUT:                                                                                                   %
%%   Xp            R x DREAMPar.d+2 - matrix with posterior samples from DREAM                              %
%%   method        string with name of marginal likelihood estimator                                        %
%%   DREAMPar      DREAM structure with algorithmic settings of DREAM                                       %
%%   Func_name     DREAM name (string between quotes) of the target distribution                            %
%%   GAME_options  (optional) GAME structure with additional options                                        %
%%     .metric     metric for optimal mixture selection (''bic''/''var'')   (default: bic)                  %
%%     .J          maximum # components mixture distribution                (def: 5)                        %
%%     .N          number of importance samples to be used                  (def: 10000)                    %
%%     .M          number of times we repeat each estimator                 (def: 1)                        %
%%     .steps      number of steps with gb or ob                            (def: 10)                       %
%%   Par_info      (optional) DREAM structure with parameter information    (default -inf/min and inf/max)  %
%%   Meas_info     (optional) DREAM structure with Meas_info data           (default = [])                  %
%%   options       (optional) DREAM structure with additional options       (default = struct)              %
%%   plug          (optional) DREAM additional input target function        (default: [])                   %
%%                                                                                                          %
%% DEFAULT:                                                                                                 %
%%   GAME_options = struct('metric','BIC','J',5,'N',1e4,'M',1,'steps',10);                                  %
%%                                                                                                          %
%% OUTPUT:                                                                                                  %
%%   Z = marginal likelihood   ( = integral of posterior pdf )                                              %
%%                             ( = integral of p(x|Y)        )                                              %
%%                             ( = integral of p(x)L(x|Y)    )                                              %
%%                                                                                                          %
%% REFERENCE:                                                                                               %
%%   E. Volpi, G. Schoups, G. Firmani, and J.A. Vrugt (2017), Sworn testimony of the model evidence:        %
%%      Gaussian Mixture Importance (GAME) sampling, Water Resources Research, 53, 6133-6158,               %
%%      doi:10.1002/2016WR020167.                                                                           %
%%                                                                                                          %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%                                                                                                          %
%% Written by Jasper A. Vrugt; with input from Elena Volpi & Gerrit Schoups                                 %
%%                                                                                                          %
%% Version 1:    June 2012       Initial setup and definition                                               %
%% Version 1.1:  Jan. 2015       Major overhaul of code                                                     %
%% Version 1.2:  Oct. 2017       Final update                                                               %
%%                                                                                                          %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%                                                                                                          %
%% MAIN IDEA:                                                                                               %
%%   Calculate the marginal likelihood as weighted mean of the ratio of the samples' target density and     %
%%   their importance density. This works because we know that importance distribution (= mix) integrates   %
%%   to 1, thus Z = 1/N*sum_i=1:N ( p(x_i)/q(x_i) ) or as we solve here in log formulation (numerically     %
%%   more stable).                                                                                          %
%%                                                                                                          %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%                                                                                                          %
%% METHODS IMPLEMENTED:                                                                                     %
%%   'ris' : reciprocal importance sampling                                                                 %
%%   'is'  : importance sampling                                                                            %
%%   'ob'  : optimal bridge sampling                                                                        %
%%   'gb'  : geometric bridge sampling                                                                      %
%%                                                                                                          %
%% NOTE 1:                                                                                                  %
%%  For geometric bridge sampling code uses: t = linspace(1,0,steps);                                       %
%%  For optimal bridge sampling code uses: nn1 = size(logQ_Xp(R1),1) * linspace(1,0.001,steps);             %
%%  If t = 0; --> RIS whereas if t = 1 --> IS; t in [0,1] bridge between both end members                   %
%%                                                                                                          %
%% NOTE 2:                                                                                                  %
%%  Make sure that you use a proper likelihood function, that is, with normalization constant included!!    %
%%                                                                                                          %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%                                                                                                          %
%% EXAMPLE:                                                                                                 %
%%   Step 1: Run example 2 of DREAM toolbox                                                                 %
%%   Step 2: Then after completion do: P = genparset(chain);                                                %
%%   Step 3: Get posterior samples: Xp = P(end-25000:end,1:DREAMPar.d+2);                                   %
%%   Step 4: Calculate marginal likelihood via GAME_sampling                                                %
%%           Z = GAME_sampling(Xp,method,DREAMPar,Func_name);                                               %
%%           where method = any option from {'ris','is','ob','gb'}                                          %
%%                                                                                                          %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
