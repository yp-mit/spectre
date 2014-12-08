function lm = feinstein_approx(n, epsil, P)
%
% 	This function computes ** APPROXIMATE ** value of Feinstein's achievability
%	bound. The normal approximation is used instead of a true distribution.
%


% conversion A->P. Old versions are all in terms of ``amplitude'' A.
A = sqrt(P);

ceps = epsil * chi2cdf(n, n);

quants = linspace(0,1,100)*ceps; quants = quants(2:end-1);

lms = [];

for qt = quants
	lgamma = inf_spect_quantile(n, qt, A);
	lm = lgamma + log2(ceps - qt);
	lms = [lms lm];

end

[lm idx] = max(lms);

%disp(sprintf(' --------> feinstein_approx(n=%d, epsil = %g, A=%g): qt_best = %.4g * ceps', ...
%	n, epsil, A, quants(idx)/ceps));

%
%	This function returns the q-th quantile of the infrmation spectrum.
%	P [ i(X^n; Y^n) <= \log gamma ] = q
%
%	TODO: we use normal approx here.
% 
function lgamma = inf_spect_quantile (n, q, A)

pp = norminv(q, 0, 1);
lgamma = pp * A/sqrt(1+A^2) * log2(exp(1)) * sqrt (n) + n/2 * log2(1 + A^2);
