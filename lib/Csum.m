function out = Csum(n,k, UBorLB, accuracy)
%sum of binomial coefficients from n choose 0 to n choose k. 
%UBorLB = 'ub' or 'lb', upper or lower bound
%accuracy: 1 - more accurate, 0 - faster

%
%   Created in 2013 by Victoria Kostina (vkostina@caltech.edu)
%

if accuracy
    out = 0;
    for j = 0:k
        out = out + C(n, j, UBorLB);
    end
else
    mult = 1;
    if strcmpi(UBorLB, 'ub')
        if n <= 2*k
            mult = Inf;
        else
            mult = (n - k)/(n - 2*k);
        end
    end
    out = mult*C(n, k, UBorLB);
end

