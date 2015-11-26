function out = Peub(R, n, d, q, g, p)
%upper bound on P_excess for coin flip source
%R - rate
%n - block length (scalar)
%d - distortion
%q - test channel parameter
%g(amma) - optimization parameter
%p - source parameter


%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

if p == 1/2
    %binary equiprobable source
    arg = (R*n - g - n*log2(2 - 2*q))/log2(q/(1-q));
    arg = max( real(arg), 0);
    %assuming d(elta) < 1/2
    out = 1 - binocdf(n*d, n, q) + binocdf( arg, n, q ) - binopdf( arg, n, q ) + exp(-2^g);
elseif p ~=1/2 %&& n <= 100
    %arbitrary binary source - exact
    %Works correctly but VERY SLOW!
    bigT = n;
    Pi = 0;
    for t = 0:bigT
        arg = (R*n - g - log2((1 - q)/(1- p))*n + t*log2((1 - q)/q))/log2((1-p)/p);
        Pi = Pi + binopdf(t,n,q)*binocdf(arg, n, p);
    end
    
    out = 2 - Pi - binocdf(n*d, n, q) + exp(-2^g);
else
    %approximate by normals
    mi = q*log2(q/(1-q)) + p*log2((1-p)/p);
    si = 1/n*(q*(1-q)*(log2(q/(1-q)))^2 + p*(1-p)*(log2((1-p)/p))^2);
    out = 2 - normcdf(R - g/n - log2((1 - q)/(1- p)), mi, sqrt(si)) - normcdf(d, q, sqrt(q*(1-q)/n)) + exp(-2^g);
end


end

