function out = binopdfbound(k, n, p, UBorLB)
%computes a fast bound to the probability of type k
%k - type number
%n - blocklength
%p - bias
%UBorLB = 'ub' or 'lb', upper or lower bound

%
%   Created in 2013 by Victoria Kostina (vkostina@caltech.edu)
% 

out = 2^(logC(n, k, UBorLB) + k*log2(p) + (n - k)*log2(1 - p));
end
