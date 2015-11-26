function out = RateDistortion(alpha, d)
%rate distortion function for binary erased source
%alpha - erasure rate of erasure channel
% d >= alpha/2 - tolerated distortion

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

p = d - alpha/2;
if p == 0
    out = 1 - alpha;
else
    out = h([(1-alpha)/2, (1-alpha)/2, alpha]) - h([1 - alpha - p, p, alpha]);
end
end