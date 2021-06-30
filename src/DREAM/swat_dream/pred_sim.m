function [ysim, ypred, yobs] = pred_sim(x, Extra)

n_model =  Extra.n_model_params;

modelparams = x(1:n_model);
errorparams = x(n_model + 1:end);

nerror = length(errorparams);

if nerror > 1
    Extra.psi = errorparams(2);
end

if nerror > 2
    Extra.lambda = errorparams(3);
end

[yobs, ysim] = run_swat(modelparams, Extra);

ysim = ysim';
yobs = yobs';

zsim = boxcoxtrans(ysim, Extra.lambda);
sigma_res = errorparams(1);
psi = Extra.psi;

Nt = length(ysim);
pd = makedist('Normal','mu',0,'sigma',sigma_res);
Wt = random(pd,1,Nt);

res = 0;
ypred = NaN(1, Nt);
for i=1:Nt
    res = psi*res + Wt(i);
    zpred = zsim(i) + res;
    ypred(i) = boxcoxinv(zpred, Extra.lambda); 
    
    if ypred(i) < 0
        ypred(i) = 0;
        res = boxcoxtrans(0, Extra.lambda) - zsim(i);
    end
end