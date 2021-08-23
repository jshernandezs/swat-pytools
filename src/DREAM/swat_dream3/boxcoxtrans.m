function z = boxcoxtrans(q, lambda)

if lambda == 0
    z = log(q);
else
    z = (q.^lambda-1)/lambda;
end