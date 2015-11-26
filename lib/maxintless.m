function x = maxintless( handle, x0 )
%maximum integer such that @handle(x) <= 0
%x0 - starting point

%
%   Created in 2013 by Victoria Kostina (vkostina@caltech.edu)
%


x = floor(x0);
while handle(x) < 0
    x = x + 1;
end
while handle(x) > 0
    x = x - 1;
end

end

