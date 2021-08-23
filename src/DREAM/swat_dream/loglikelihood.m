function ret = loglikelihood(x, Extra)

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

zobs = boxcoxtrans(yobs, Extra.lambda);
zsim = boxcoxtrans(ysim, Extra.lambda);

res = zobs - zsim;

Nt = length(res);
sigma_res = errorparams(1);
psi = Extra.psi;

pd = makedist('Normal','mu',0,'sigma',sigma_res);

ret = -Nt * log(sigma_res) - Nt/2 * log(2*pi) ...
      - 1/(2*sigma_res^2) * sum((res(2:end) - psi * res(1:end-1)).^2)...
      + sum(log(boxcoxJ(yobs, Extra.lambda)))...
      - sum(log(1 - cdf(pd, boxcoxtrans(0, Extra.lambda) - zsim(2:end) - psi * res(1:end-1)) + realmin));
