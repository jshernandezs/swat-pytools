clc, clear, close all

% settings
data = load('po_prior.csv');
bw = [0.0065,0.0066,1.89,0.017,0.018,0.27,0.021,0.83,0.0041,1.23,0.017,0.00097,0.059,0.030,9.7e-05,0.0043];

X = logmvkrnd(data,bw,100000);

f = mvkpdf(X,data,bw);
logf = log(f);

logf2 = logmvkpdf(X,data,bw);