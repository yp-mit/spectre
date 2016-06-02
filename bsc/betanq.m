function lb = betanq(nq,n, delta, shutup)
% returns log \beta_(1-nq) for BSC(\delta)^n
%
% Note we pass nq instead of q to save precision as q = 1-eps+tau usually

if (nargin < 4) || (max(size(shutup)) == 0)
	shutup = 0;
end

% Step 1. Find correct gamma in Neyman-Pearson

if(delta > .5)
	disp('Does not work for delta > .5');
	error('betanq');
end

%
% Overshoot = P_{Y|X} (dP/dP >= gamma) - q
%
if( nq == 1 )
	T = 0;
	overshoot = binopdf(0,n,delta);
else
	T = n;
	ncumdist = 0;
	overshoot = nq;

	while 1
		ncumdist = ncumdist + binopdf(T, n, delta);
		if ncumdist > nq
			break;
		end
		T = T - 1;
		overshoot = nq - ncumdist;
	end
end

if ~shutup
	disp(sprintf('------ betanq(nq = %g, n = %d, delta = %.2f) -------', nq, n, delta));
	disp(sprintf('     Step 1 : T = %d, T/n = %g', T, T/n));
end

log_gamma = log2(delta/(1-delta)) * T + n * log2(2*(1-delta));

if ~shutup
	disp(sprintf('            : overshoot = %g, log_gamma = %g', overshoot, log_gamma));
end

% Paranoia
if(overshoot < 0), disp('error: overshoot < 0'); error('betanq'); end

%
% Step 2. Compute term1 and term2/term1
%
% Term1 = P_Y{dP/dP <= \gamma} = sum_{k=0}^T nchoosek(n,k) 2^{-n}
%
% Term2 = -overshoot / gamma;
%
%
% We compute log(Term1) using our log-sum method. Precision of this
% method is approximately 2^(-1024) * (number of summands). So it is
% much less than eps = 2^(-52).
%
% Indeed, we compute (assume A>B)
%      log(A + B) = log A  + log (1 + B/A);
% B/A = 2^{log B - log A}. So we only lose precision when B/A gets below 2^{-1024}.
% But log (1 + x) <= x. Thus each time we get error less than 2^{-1024}.
%

log_sum = -Inf;
for k = 0:T;
	cur_term = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1);
	if cur_term > -Inf
		if cur_term > log_sum
			log_sum = cur_term + log1p(exp(log_sum - cur_term));
		else
			log_sum = log_sum + log1p(exp(cur_term - log_sum));
		end;
	end
end
% correct for the log base and correct for 2^{-n} term
log_sum = log_sum / log(2) - n;

% term1 = 2^log_sum

log_ratio = log2(overshoot) - log_sum - log_gamma;

% More paranoia
if(log_ratio > 5*eps)
	disp('Error: Can not happen: log_ratio < 0');
	error('betanq');
end

if  (log_ratio > 0) && (log_ratio <= 5*eps) 
	log_ratio = 0;
end


if ~shutup
	disp(sprintf('     Step 2 : ratio = %g      ', 2^log_ratio));
end

%
% Step 3. Compute beta_q
%

lb = log_sum + log1p(-2^log_ratio)/log(2);

if ~shutup
	disp(sprintf('     Step 3 : lb = %g      ', lb));
end
