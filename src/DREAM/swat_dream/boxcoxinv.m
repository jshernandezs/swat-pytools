function q = boxcoxinv(z, lambda)

if lambda == 0
    q = exp(z);
else
    q = (lambda * z + 1).^(1/lambda);
end