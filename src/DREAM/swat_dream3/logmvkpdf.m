function logf = logmvkpdf(x, data, bw)

[n,d] = size(data);

s = 0;
for i=1:n
    p = 0;
    for j=1:d
        p = ((x(:,j) - data(i,j)) / bw(j)) .^ 2  + p;
    end
    s = s + exp(-0.5*p);
end

logf = -log(n*(2*pi)^(d/2)*prod(bw)) + log(s);

if isinf(logf)
    logf = -(realmax);
end