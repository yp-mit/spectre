function lm = wolfowitz(n, epsil, P)
% returns wolfowitz' upper bound on Log M

% conversion A->P. Old versions are all in terms of ``amplitude'' A.
A = sqrt(P);
q = 1 - epsil;


%
% Wolfowitz bound
%

deltas = (1-q).* [ .01 .05 linspace(.1, .3, 10) .4 .5 .6 .7 .8 .9 1];
%deltas = (1-q) .* linspace(0,1,100); deltas = deltas(2:end-1);
Qs = q - deltas;
lbeta_last = -Inf;
delta_last = [];

%
% TODO: note that if we parameterize over pp0 instead of qc we could eliminate
%	very costly ncx2inv() from the loop
%

for qc = Qs;
		pp0 = ncx2inv(qc, n, n/A^2); 
		delta = q - ncx2cdf(pp0, n, n/A^2); 
		gammatil = (1 + A^2) * n - A^2 * pp0; 
		lgamma = gammatil * log2(exp(1)) / (2  + 2*A^2) + n/2 *log2 (1 + A^2); 
		if(delta < 0)
			disp('ERROR: delta should never be zero');
			error('betaq_low_v2');
		end
		lbeta_est = log2(delta) - lgamma;
		if(lbeta_est < lbeta_last)
			break;
		end

		delta_last = delta;
		lbeta_last = lbeta_est;
end

    disp(sprintf(['wolfowitz(n=%d, epsil = %g, P=%g): WOLFOWITZ: '...
	'delta_best = %.2f, log M = %.1f'], ...
	n, epsil, P, delta_last / (1-q), -lbeta_last));

lm = -lbeta_last;
