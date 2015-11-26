function out = loggamma(x, UBorLB)
%natural logarithm of the gamma function
%UBorLB = 'ub' or 'lb', upper or lower bound

%
%   Created in 2013 by Victoria Kostina (vkostina@caltech.edu)
%


out = gamma(x);
if isfinite(out)
    out = log(out);
    return;
end
%  Sharp upper and lower bounds for the gamma function
% Horst Alzer
beta = 1/1620;
y = x - 1;
switch lower(UBorLB)
    case 'ub'
        out = .5*log(2*pi*y) + y.*log(y/exp(1)) + y/2.*log(y.*sinh(1/y)).*(1 + beta/y^5); %upper bound
    case 'lb'
        out = .5*log(2*pi*x) + y.*log(y/exp(1)) + y/2.*log(y.*sinh(1/y)); %lower bound
    otherwise
        disp('gammafun(x, UBorLB): UB or LB?');
        out = NaN;
end
end

