function f = mdist(n, delta, Gammas)
% returns the log2 of value of P[i(X,Y) <= \gamma] + 2^\gamma P[i(X,\bar Y) > \gamma] for BEC

f = zeros(1,size(Gammas,2));
Ks = 0:n;
terms0 = (gammaln(n+1) - gammaln(Ks+1) - gammaln(n-Ks+1))/log(2) + n*log2(delta) + Ks*log2((1-delta)/delta);

for idx = 1:size(Gammas,2);
	gamma = Gammas(idx);

	terms = terms0 + min(gamma - Ks, 0);

	f(idx) = sumlog2(terms);
end
