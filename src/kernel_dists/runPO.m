% Sampling model results from Pareto-optimal results
clc, clear, close all
rng('default')

%% settings
inputFolder = 'Data';
filename = 'filtered_nds_params.csv';
filename2 = 'po_simulations.csv';
nprocessors = 4;

plugin.n_model_params = 15;
plugin.path_to_python = '~/anaconda3/envs/swatpy/bin/python';
plugin.lambda = 0.2;                                                            % Box-cox transformation parameter
plugin.psi = 0.8;                                                               % Autocorrelation error model parameter
plugin.wrapper = 'swat_wrapper.py';

addpath('../DREAM/swat_dream')
%% get parameter sets
d = 16;
n_model_param = plugin.n_model_params;

ParEval = load([inputFolder '/' filename]);
Neval = length(ParEval);

ParEval = [ParEval ones(Neval,1)];
%% run predictions (parameter and total error)

% run the first model (we need it to get number of time steps)
[model, ~, obs] = pred_sim(ParEval(1, :), plugin);
Nt = length(obs);

% models = zeros(Neval, Nt);
% models(1,:) = model;
% 
% % execute model runs
% parpool(nprocessors);
% 
% parfor i=2:Neval
%     [models(i,:), ~] = pred_sim(ParEval(i, :), plugin);
% end
% 
% delete(gcp('nocreate'))

models = load([inputFolder '/' filename2]);

%% generate distribution for sigma_w
sigma_w = NaN(Neval,1);
for i=1:Neval
    sim = models(i,:);
    res = boxcoxtrans(obs, plugin.lambda) - boxcoxtrans(sim, plugin.lambda);
    wt = res(2:end) - plugin.psi*res(1:end-1);
    sigma_w(i) = std(wt);
end

save('sigma_w2','sigma_w')

%% generate predictions

preds = models;
psi = plugin.psi;

for j=1:Neval
    ysim = models(j,:);
    zsim = boxcoxtrans(ysim, plugin.lambda);
    
    Nt = length(ysim);
    pd = makedist('Normal','mu',0,'sigma',sigma_w(j));
    Wt = random(pd,1,Nt);

    pd0 = makedist('Normal','mu',0,'sigma',sigma_w(j)/sqrt(1 - psi^2));
    res = random(pd0);

    ypred = NaN(1, Nt);
    for i=1:Nt
        if i > 1
            res = psi*res + Wt(i);
        end
        zpred = zsim(i) + res;
        ypred(i) = boxcoxinv(zpred, plugin.lambda); 

        if ypred(i) < 0
            ypred(i) = 0;
            res = boxcoxtrans(0, plugin.lambda) - zsim(i);
        end
    end

    preds(j,:) = ypred;
end

save('PO_simulations2','models','preds')
%% plots and performance

PredInt = 95; % prediction uncertainty range

[N_Pars, ~] = size(preds); % number of samples 

Lb = floor( (100 - PredInt)/200 * N_Pars); % define lower bound of prediction interval
Ub = N_Pars - Lb; % define upper bound of prediction interval

% initialize output variables
Fq = NaN(1,Nt);
cv = Fq;
par_unc = NaN(Nt,2);
tot_unc = par_unc;

% obtain quantiles per time step

for i=1:Nt
    y = preds(:, i);
    [f, x] = ecdf(y); % empirical cdf for ith time step
    Fq(i) = interp1(x(2:end),f(2:end),obs(i),'linear',0); % quantile for ith observation
    cv(i) = std(y)/mean(y); % coef. of variation for ith time step
    
    % Sort model output from low to high
    ymod = models(:, i);
    a = sort(ymod);
    % And take the desired prediction uncertainty ranges
    par_unc(i,:) = [a(Lb), a(Ub)]; 
    % Same with total uncertainty
    a = sort(y);
    % And take the desired prediction uncertainty ranges
    tot_unc(i,:) = [a(Lb), a(Ub)];    
end

pd = makedist('uniform'); % create uniform [0,1] distribution
[pqq, u] = ecdf(Fq); % obtain empirical cdf from observed quantiles

Fu = cdf(pd, Fq);
Fs = interp1(u(2:end),pqq(2:end),Fq,'linear',0);

Rqq = 2/Nt * (sum(abs(Fu - Fs)));
Pqq = mean(cv);

ypred_mean = mean(preds, 1);
Vqq = abs((sum(obs) - sum(ypred_mean))/sum(obs));

figure()
hold on
plot([0 1], [0 1], 'k--')
plot(u, pqq)
axis square
box on
xlabel('Quantile U[0,1]')
ylabel('Obs. quantile')

figure('units','normalized','position',[0.25 0.25 0.4 0.4])
hold on
Fill_Ranges((0:Nt-1),tot_unc(:,1),tot_unc(:,2),[0.85 0.85 0.85]);
Fill_Ranges((0:Nt-1),par_unc(:,1),par_unc(:,2),[0.4 0.4 0.4]);
plot(obs, 'r--')
xlabel('Time (days)')
ylabel('Streamflow (m^3 s^{-1})')
legend('Total uncertainty', 'Parameter uncertainty', 'Observations','Location','NW')
xlim([-Inf Inf])
box on