function [lbeta lgamma] = betaq_low_v2(q, n, P)
%
% This a version 2 of betaq_low. See comment in betaq_low.m for the idea.
%
% Plan:
%	1. use -log \beta_q ~~  nC + K \sqrt(n) to estimate the value of beta. 
%		NOTE: Trials in script_optimize_wolfowitz.m showed that
%		in fact ncx2cdf() starts returning zero at around 2^-570 => 
%		we can rely on the output of ncx2cdf() to the point when
%		it returns _ZERO_.
%
%		So a simple test (term1 == 0) is OK for switching b/w method 2 and 3:
%`	2. If it is much higher than minimum double exponent -- use old method: log2(term1 + term2)
%	3. Otherwise use the wolfowitz bound optimization.
%
% Corrections:
%	-- For the old method:  subtract 2*eps*q from delta before computing lower bound
%	-- For the optimization of the Wolfowitz:
%		should aim at value delta ~ 0.1..0.3 (1-q)
%		I.e. try setting Fp \approx q - 0.2*(1-q)
%

%
%	Fp(gamma) = P[ dP/dQ >= gamma ] = P(  sum (Z_i - 1/A)^2 <= pp  )
%	Fq(gamma) = Q[ dP/dQ >= gamma ] = P(  sum (Z_i - sqrt(1+A^2)/A)^2 <= qq )
%
% so, Fp(2^lgamma) = ncx2cdf(pp0, n, n/A^2);
%     Fq(2^lgamma) = ncx2cdf(qq0, n, n*(1+1/A^2));
%
%   and this pair (Fp,Fq) gives us a point on the ROC curve (q', beta_q')
%   and also we tried very hard to choose gamma parameter so that this particular point
%   (Fp, Fq) is very-very close to (q, beta_q) (remember we are seeking beta_q at a particular 
%   point).
%
%   However, since Fp != q  just to be on the safe side we are using a lower bound on beta_q
%   that is constructed from any point (Fp, Fq) and knowing that the slope at this point 1/gamma.


% conversion A->P. Old versions are all in terms of ``amplitude'' A.
A = sqrt(P);

%
% First we try to compute betaq precisely (almost)
%
pp0 = ncx2inv(q, n, n/A^2);
delta = q - ncx2cdf(pp0, n, n/A^2);
gammatil = (1 + A^2) * n - A^2 * pp0;
lgamma = gammatil * log2(exp(1)) / (2  + 2*A^2) + n/2 *log2 (1 + A^2);
qq0 = ((1+A^2) * n - gammatil) / ((1+A^2)*A^2);

%
% for pp0 that is very close to true q-th quantile of P[ dP/dQ >= gamma ] we
% will have that term1 is very very close to true beta. 
%
% If term1 = 0 -- then use another method
%
term1 = ncx2cdf(qq0, n, n*(1 + 1/A^2));
term2 = (delta - 2*q*eps) * 2^(-lgamma);

if (term1 > 0)
	lbeta_prec = log2(term1+term2);

	disp(sprintf(['betaq_low_v2(q=%g, n=%d, P=%g): PRECISE: '...
		'term2/term1 = %.1g, delta=%.1g, lgamma = %.3g, lbeta = %.1f'], ...
		q,n,P, abs(term2/term1), delta, lgamma, lbeta_prec));

	lbeta = lbeta_prec;
else

	%
	% Compute using new ncx2log()
	%
	logterm1 = ncx2log(qq0, n, n*(1+1/A^2)) * log2(exp(1));
    if( delta > 2*q*eps)
    	logterm2 = log2(delta - 2*q*eps) - lgamma;
    else
        logterm2 = -Inf;
    end

	if(logterm1 > -1240) 

		if (logterm1 > logterm2)
			lbeta = logterm1 + log2(1+2^(logterm2 - logterm1));
		else
			disp('should not happen, investigate!');
			error('betaq_low_v2');
		end

		disp(sprintf(['betaq_low_v2(q=%g, n=%d, P=%g): LOG_PRECISE: '...
			'term2/term1 = %.1g, delta=%.1g, lgamma = %.3g, lbeta = %.1f'], ...
		q,n,P, 2^(logterm2 - logterm1), delta, lgamma, lbeta));
		return;
	end

	%
	% Wolfowitz approximation 
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
	disp(sprintf(['betaq_low_v2(q=%g, n=%d, P=%g): WOLFOWITZ: '...
			'delta_best = %.2f, lbeta = %.1f'], ...
		q,n,P, delta_last / (1-q), lbeta_last));
	lbeta = lbeta_last;
end
