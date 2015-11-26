function out = PeubShannon(R, n, d, q, g)
%upper bound on P_excess for Gaussian source
%R - rate
%n - block length (scalar)
%d - distortion
%q - test channel parameter
%g(amma) - optimization parameter

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

BIG = 3*4; %3*variance
arg = (n*R - g - n/2*log2(1/(q^2)))/log2(exp(1))*2;
if n <= 100
out = 2 - diffChi2Cdf( arg, n) - chi2cdf(n*d/(q^2),n) + exp(-2^g);
else
    %approximate by Gaussians:
    out = 2 - normcdf( arg/sqrt(n)/2, 0,1) - chi2cdf(n*d/(q^2),n) + exp(-2^g);
end


    function out = diffChi2Cdf(a, n)
         out = n*quad(@density, 0, BIG);
         function out = density(r)
             out = chi2pdf(r*n, n).*chi2cdf(max(0, r*n + a), n);
         end
    end        
end

