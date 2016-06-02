function [lm lm_avg] = rcu_ach(n, delta, epsil, plow, pup);
% Compute RCU achievability bound. 
% lm returns the bound valid for maximal probability of error
% lm_avg is for the average
%
% This is not a vectorized version, because the function is really slow.
% You can speed it up significantly if you provide a bracket plow, pup for the value of logm.
% For example, a smart idea is to set 
% plow = gallager_ach(n, delta, epsil) and pup = converse(n,delta,epsil);

if nargin < 5
	plow = 0;
	pup = n;
end

eps_test = precise_rand(n, delta, pup);
if(eps_test < epsil)
	disp(sprintf([	'-- achiev_prand(n = %d, delta = %g, epsil = %g): eps_test = %g\n'...
			'        This is a bug? precise_rand() contradicts converse ?!?!'], ...
			n, delta, epsil, eps_test));
	error('achiev_prand');
	lm_avg = 0; lm = 0;
	return;
end

eps_test = precise_rand(n, delta, plow);
if(eps_test > epsil)
	disp(sprintf([	'-- achiev_prand(n = %d, delta = %g, epsil = %g): eps_test = %g\n'...
			'        Can not find lower bound for precise_rand() !!!!'], ...
			n, delta, epsil, eps_test));
	%error('achiev_prand');
	lm_avg = 0; lm = 0;
	return;
end

% Take into account that we are computing AVERAGE prob. of error,
% and so we need a random linear code trick to go to maximum => log M must be integer
while floor(plow) < floor(pup);
	ptest = (plow + pup)/2;
	eps_test = precise_rand(n, delta, ptest);
	if(eps_test > epsil)
		pup = ptest;
	else
		plow = ptest;
		% This is not needed as we stop on floor-test
		%if  (epsil-eps_test)< 1e-2*epsil
		%	break;
		%end
	end
end
lm_avg = plow;
lm = floor(plow);
