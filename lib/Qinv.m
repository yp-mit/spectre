function out = Qinv(x) 
%inverse Q-function 
%
%   Created in 2013 by Victoria Kostina (vkostina@caltech.edu)
%
        out = sqrt(2)*erfcinv(2*x);
end