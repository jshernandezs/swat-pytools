% Get performance indicators and plots
clc, clear, close all

filename = 'MTDREAM_Results/predictions.mat';
load(filename)

PredInt = 95; % prediction uncertainty range

Nt = length(obs); % number of observations
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
ylabel('Streamflow ($m^3 s^{-1}$)')
legend('Total uncertainty', 'Parameter uncertainty', 'Observations','Location','NW')
xlim([-Inf Inf])
box on