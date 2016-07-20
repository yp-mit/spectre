function logm=converse_spec2(n, delta, epsil, lm_min, lm_max, shutup)
% This is a special converse for BEC, based on channel-realization analysis. 
% This is a non-vectorized version!
% You may provide a bracket lm_min, lm_max to speed up the search.
% This bound is for average P_e.

if(nargin < 6) || isempty(shutup)
	shutup = 0;
end

if (nargin < 5) || (max(size(lm_min))==0) || (max(size(lm_max))==0)
	%% First approximation
	%%	we hope that this is a lower bound on log M
	lm_min = binoinv(epsil, n, 1-delta);
	if(epsil < .5)
		lm_max = n*(1-delta);
	else
		lm_max = n;
	end
end

if(~shutup)
	disp(sprintf('-- converse_spec2(n=%d, delta=%g, epsil=%g): lm_min = %.2f, lm_max = %2.f', ...
		n, delta, epsil, lm_min, lm_max));
end

Ls = 0:n;
bino_coeffs = (gammaln(n+1) - gammaln(Ls + 1) - gammaln(n - Ls + 1))/log(2) + Ls.*log2(delta) + (n-Ls).*log2(1-delta);

% Check lm_min and lm_max
L0 = floor(n - lm_min) + 1;
terms = bino_coeffs(L0+1:end) + log1p(-2.^(n - Ls(L0+1:end) - lm_min))/log(2);
if(2^sumlog2(terms) > epsil)
	disp('ERROR: lm_min is not a lower-bound!');
	error('lm_min');
end
L0 = floor(n - lm_max) + 1;
terms = bino_coeffs(L0+1:end) + log1p(-2.^(n - Ls(L0+1:end) - lm_max))/log(2);
if(2^sumlog2(terms) < epsil)
	disp('ERROR: lm_max is not an upper-bound!');
	error('lm_max');
end

nr_iter = 1;

while 1;
	lm_test = (lm_min + lm_max)/2;
	L0 = floor(n - lm_test) + 1;
	terms = bino_coeffs(L0+1:end) + log1p(-2.^(n - Ls(L0+1:end) - lm_test))/log(2);
	test_epsil = 2^sumlog2(terms);
	if(test_epsil < epsil)
		if(~shutup)
			disp(sprintf('    Step %d. [%.2f %.2f] lm_test = %.2f, test_epsil = %.2g, heading UP', ...
				nr_iter, lm_min, lm_max, lm_test, test_epsil));
		end
		lm_min = lm_test;
	else
		if(~shutup)
			disp(sprintf('    Step %d. [%.2f %.2f] lm_test = %.2f, test_epsil = %.2g, heading DOWN', ...
				nr_iter, lm_min, lm_max, lm_test, test_epsil));
		end
		lm_max = lm_test;
	end
	
	if (lm_max - lm_min) < 1e-2 
		break;
	end
	nr_iter = nr_iter + 1;
end


if(~shutup)
	disp(sprintf('  ---> log M <= %.2f', lm_max));
end

logm = lm_max;
