function [lbeta lgamma method] = betaq_up_v2(q, n, P)
% computes the _UPPER_ bound (as good as it can find) for the beta_q function:
%   beta_q = min Q(B), s.t. P(B) >= q
%
%   where 
%
%	P ~ Gsn( x_0,	I_n )
%	Q ~ Gsn( 0, (1+A^2) I_n )
%
%	And x_0 something that is ||x_0||^2 = n A^2
%
%	Taking x_0 = [A, ... A] we get
%
%	P[ dP/dQ >= gamma ] = P(  sum (Z_i - 1/A)^2 <= pp  )
%	Q[ dP/dQ >= gamma ] = P(  sum (Z_i - sqrt(1+A^2)/A)^2 <= qq )
%
%	pp = ((1 + A^2) n - gammatil) / A^2
%	qq = ((1 + A^2) n - gammatil) / ((1+A^2) A^2)
%
%	gammatil = [ log gamma - n/2 log (1+A^2) ] 2 (1 + A^2) / log e
%
%
% Idea for the upper bound: take any \gamma s.t. P_\gamma >= q. Then Q_\gamma >= \beta_q
% Or for high gamma: \beta_q <= 1/gamma;


% conversion A->P. Old versions are all in terms of ``amplitude'' A.
A = sqrt(P);

pp0 = ncx2inv(q, n, n/A^2);
iter = 1;
while 1; 
	pgam = ncx2cdf(pp0, n, n/A^2);

	if (pgam >= q)
		break
	end

	delta = q - pgam;
	%
	% For q > .5 this is underestimate => we will undershoot.
	%
	pp0 = pp0 + delta/ncx2pdf(pp0, n, n/A^2);
	iter = iter + 1;
end

delta = pgam - q;


gammatil = (1 + A^2) * n - A^2 * pp0;
lgamma = gammatil * log2(exp(1)) / (2  + 2*A^2) + n/2 *log2 (1 + A^2);
qq0 = ((1+A^2) * n - gammatil) / ((1+A^2)*A^2);

%
% The funny property of ncx2cdf is that it returns pretty consistent results up until 
% it starts returning 0's. So we believe it until term0. After that we use simple 1/gamma
% upper bound. 
%
% TODO: change 1/gamma bound to a large deviation chernoff-type one.
%	04/18/2007: done. (see various comparisons in script_optimize_upper)
%
% 05/02/2007: all upper bounding is deprecated, now we simply use ncx2log() which
%	works for all values of interest.
%
term1 = ncx2cdf(qq0, n, n*(1 + 1/A^2));

if (term1 == 0) 
	%
	% This is the old upper bound beta <= 1/gamma
	%
	%lbeta =  -lgamma;


	lbeta = log2(exp(1)) * ncx2log(qq0, n, n*(1 + 1/A^2));
    
    if (lbeta > -1240)
       	if (lbeta > -lgamma)
    		disp('ERROR: precise value of lbeta is above 1/gamma ');
        	error('betaq_up_v2');
    	end
        msg = 'PREC_LOG';
    else
    
    %
	% LARGE-DEVIATIONS method: if everything else fails
	%
	pstar = A^2/2; 
	sn = qq0/n;

	s = warning('off', 'optim:fminunc:SwitchingMethod');
	opts = optimset('Display', 'off');
	[popt fopt] = ...
		fminunc(@(p)( -1/2 * log(2*abs(p)+1) - abs(p)/(2*abs(p)+1) * (1+ 1/A^2) + abs(p)*sn ), ...
			pstar, opts);
	warning(s);

	lbeta_ld = fopt * log2(exp(1)) * n;
	lbeta_gamma = -lgamma;

	if (lbeta_gamma > lbeta_ld)
		lbeta = lbeta_ld;
		msg = 'LARGE-DEV';
		method = 2;
	else
		lbeta = lbeta_gamma;
		msg = 'GAMMA';
		method = 1;
	end;
    
    end;

	disp(sprintf('betaq_up(q = %g, n = %d, A = %g): %s: iter = %d, delta = %.1g, lbeta = %.2f', ...
			q, n, A, msg, iter, delta, lbeta));    
else
	lbeta = log2(term1);
	method = 0;
	disp(sprintf('betaq_up(q = %g, n = %d, A = %g): PRECISE: iter = %d, delta = %.1g, lbeta = %.2f (loose: %.1f %% worse) ', ...
			q, n, A, iter, delta, lbeta, 100*(-lgamma-lbeta)/(-lbeta)));
end;


