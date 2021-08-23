function R = mvkrnd(data, bw, n)

if nargin < 3
    n = 1;
end

gm = gmdistribution(data, bw.^2);
R = random(gm, n);