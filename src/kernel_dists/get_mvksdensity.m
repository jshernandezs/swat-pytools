% Kernel density estimation
clc, clear, close all
rng('default')

inputFolder = 'Data';
filename = 'filtered_nds_params.csv';
k = 5; % number of partitions for k-fold cross validation

data = load([inputFolder '/' filename]);
load('sigma_w.mat')

data = [data, sigma_w];

% Silverman's rule of thumb
[n, d] = size(data);

bw = NaN(1,d);
for i=1:d
    y = data(:,i);
    sdev = std(y);
    bw(i) = sdev*(4/((d+2)*n))^(1/(d+4));
end

%lb = 0.10 * bw;
%ub = 1.10 * bw;

%cv = cvpartition(n, 'Kfold', k);

%z = dist_cdf(bw, data, cv);

%options = optimoptions('ga','Display','iter','UseParallel',true,...
%                   'MaxStallGenerations', 50, 'FunctionTolerance', 1e-6);
               
%tic               
%bw_opt = ga(@(x)dist_cdf(x, data, cv), d, [], [], [], [], lb, ub, [], options);
%toc

load('bw.mat')

gm = gmdistribution(data, bw_opt.^2);
params = random(gm, 10000);

figure(1)
for i=1:d
    subplot(4,4,i)
    hold on
    histogram(data(:,i), 15, 'facecolor', 'b', 'Normalization', 'pdf')
    histogram(params(:,i), 15, 'facecolor', 'r', 'Normalization', 'pdf')
end

figure(2)
for i=1:d
    [f1, x1] = ecdf(data(:,i));
    [f2, x2] = ecdf(params(:,i));
    subplot(4,4,i)
    hold on
    plot(x1, f1, 'b')
    plot(x2, f2, 'r')
end

% save('bw','bw_opt')