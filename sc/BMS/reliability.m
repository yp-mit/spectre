function out = reliability(R, d)
%reliability function for binary source
%R - rate
%d - excess distortion

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

%source:
global gP;
p = gP;

%optimum distribution
q = real( fsolve(@(q)h(q) - h(d) - R,.11)); 

%relative entropy:
out = D([q 1-q], [p 1-p]);

end

