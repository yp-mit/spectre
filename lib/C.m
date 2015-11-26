function out = C(n,k, UBorLB)
%computes fast bounds to the binomial coefficient
%UBorLB = 'ub' or 'lb', upper or lower bound

%
%   Created in 2013 by Victoria Kostina (vkostina@caltech.edu)
%
if n < k || k < 0
    out = 0;
elseif n < 500
    out = nchoosek(n,k);
elseif k == 0 || n - k == 0
    out = 1;
else
    switch lower(UBorLB)
        case 'ub'
            out = sqrt(n/(2*pi*k*(n-k)))*2^(n*h(k/n)); %upper bound
        case 'lb'
            out = sqrt(n/(8*k*(n-k)))*2^(n*h(k/n)); %lower bound
        otherwise
            disp('Error: C(n,k, UBorLB): UB or LB?');
            out = NaN;
    end
    %out = 2^(n*h(k/n));
end
end

