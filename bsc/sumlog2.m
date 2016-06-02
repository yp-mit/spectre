function [f cumf] = sumlog2(S)
% returns log2(sum(2.^S)) in a stable way

% We compute log2(f) using our log-sum method. Precision of this
% method is approximately 2^(-1024) * (number of summands). So it is
% much less than eps = 2^(-52).
%
% Indeed, we compute (assume A>B)
%      log(A + B) = log A  + log (1 + B/A);
% B/A = 2^{log B - log A}. So we only lose precision when B/A gets below 2^{-1024}.
% But log (1 + x) <= x. Thus each time we get error less than 2^{-1024}.
%

const_log2e = 1/log(2);
log_sum = -Inf;
cumf = zeros(size(S));
idx = 1;
for cur_term = S;
	if cur_term > -Inf
		if cur_term > log_sum
			log_sum = cur_term + log1p(2^(log_sum - cur_term))*const_log2e;
		else
			log_sum = log_sum + log1p(2^(cur_term - log_sum))*const_log2e;
		end;
	end
	cumf(idx) = log_sum;
	idx = idx + 1;
end

f = log_sum;
