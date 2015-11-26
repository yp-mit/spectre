function out = PelbMarton(R, n, d, q, delta, p)
%lower bound on P_excess for coin flip source
%R - rate
%n - block length (scalar)
%d - excess distortion
%q - auxiliary source distribution
%delta - "radius of looseness"
%p - source distribution parameter


%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

%relative entropy:
Dav = D([q 1-q], [p 1-p]);

%distortion-rate function
DR = real( fsolve(@(x)h(q) - h(x) - R,.11)); 

%constants:
a = (DR - d)/(1 - d);
A = delta/(log2(q*(1-p)/(p*(1-q))));

out = (a + binocdf((q + A)*n, n, q)  - 1 )*(2)^(-n*(Dav + delta));


end

