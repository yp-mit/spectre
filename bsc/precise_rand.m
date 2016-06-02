function epsil = precise_rand(n, delta, logm, shutup)
% This function computes (almost) precise random coding bound for BSC. The only imprecision is due to 
% union bound.
% epsil <= sum_k  \delta^k (1-delta)^(n-k) P[one of M-1 other codewords is inside sphere <= k]

Ks = 0:n;
bino_coeffs = (gammaln(n+1) - gammaln(Ks+1) - gammaln(n-Ks+1))/log(2);
terms_spheres = bino_coeffs + logm - n;
[tmp log_spheres] = sumlog2(terms_spheres);

% cut off stupid >1 values (this what gallager rho-trick is supposed to do also)
log_spheres = min(log_spheres, 0);

terms = bino_coeffs + log2(delta)*Ks + log2(1-delta)*(n-Ks) + log_spheres;

epsil = 2^sumlog2(terms);

if(nargin < 4) || (isempty(shutup))
	disp(sprintf('-- precise_rand(n = %d, delta = %g, logm = %.1f): epsil = %.3g', n, delta, logm, epsil));
end
