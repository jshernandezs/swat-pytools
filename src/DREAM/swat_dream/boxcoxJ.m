function z = boxcoxJ(q, lambda)
if lambda == 0
    z = 1./q;
else
    z = q.^(lambda-1);
end