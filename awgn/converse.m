function Lms = converse(Ns, epsil, P)
% Compute meta-converse upper bound on \log M^*(n, \epsilon, P)
Lms = [];
for n = Ns;
	Lms = [Lms -betaq_low_v2(1-epsil, n+1, P)];
end
