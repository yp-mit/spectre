function out = logC(n,k, UBorLB)
%binary logarithm of the binomial coefficient
%binary log of C(n,k)
%UBorLB = 'ub' or 'lb', upper or lower bound

%
%   Created in 2013 by Victoria Kostina (vkostina@caltech.edu)
%
if n < k || k < 0
    out = -Inf;
elseif n < 500
    out = log2(nchoosek(n,k));
elseif k == 0 || n - k == 0
    out = 0;
else
    switch lower(UBorLB)
        case 'ub'
            out = log2(sqrt(n/(2*pi*k*(n-k)))) + n*h(k/n); %upper bound
        case 'lb'
            out = log2(sqrt(n/(8*k*(n-k)))) + n*h(k/n); %lower bound
        otherwise
            disp('Error: C(n,k, UBorLB): UB or LB?');
            out = NaN;
    end
end
end

