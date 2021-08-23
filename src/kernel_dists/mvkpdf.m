function f = mvkpdf(x, data, bw)

gm = gmdistribution(data, bw.^2);
f = pdf(gm, x);