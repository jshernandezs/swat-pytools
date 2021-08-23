function ret = dist_cdf(u, obs, cv)

[~, d] = size(obs);
z = NaN(1,cv.NumTestSets);

for i=1:cv.NumTestSets
    
    training = obs(cv.training(i),:);
    test = obs(cv.test(i),:);
    gm = gmdistribution(training, u.^2);
    sim = random(gm, 10000);
    
    dd = NaN(1,d);
    for j=1:d
        [fobs, xobs] = ecdf(test(:,j));
        [f, x] = ecdf(sim(:,j));

        fsim = interp1(x(2:end),f(2:end),xobs,'pcubic',0);
        dd(j) = sum((fobs - fsim).^2);    
    end

    z(i) = mean(dd);    
end

ret = mean(z);