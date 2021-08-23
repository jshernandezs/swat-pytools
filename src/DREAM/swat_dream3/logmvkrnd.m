function R = logmvkrnd(data, bw, n)

if nargin < 3
    n = 1;
end

gm = gmdistribution(data, bw.^2);
R = random(gm, n);

% truncate values

lb = [0,-0.25,0,0,0,0,0,0,0.02,0,0,0,0,-0.25,1,0.1];
lu = [1,0.25,100,1,1,500,1,5000,0.2,1000,1,0.3,500,0.25,24,5];

[Nt,d] = size(R);

for i=1:d
    for j=1:Nt
        if R(j,i) < lb(i)
            R(j, i) = lb(i) + abs(lb(i) - R(j,i));
        end
        
        if R(j,i) > lu(i)
            R(j, i) = lu(i) - abs(lu(i) - R(j,i));
        end
    end    
end