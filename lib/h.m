function out = h(p)
%entropy (base 2)

%
%   Created in 2013 by Victoria Kostina (vkostina@caltech.edu)
%

if p == 0
    out = 0;
    return;
end

s = sum(p);
if s < 1
    p = [p 1-s];
end
out =  -sum(p.*log2(p));

end

