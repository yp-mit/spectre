function [out fval]= fzerointeger(fhandle, startk0, endk, type)
%finds max/min k such that fhandle(k) <= 0
%type = 'min' or 'max'

%
%   Created in 2013 by Victoria Kostina (vkostina@caltech.edu)
%



fcache = NaN(1, endk - startk0 + 1);
startk = startk0;
stepk = ceil((endk - startk + 1)/8);

while 1
    [out fval] = bestongrid(stepk, startk, endk);
    
    switch lower(type)
        case 'min'
            startk = max(out - stepk, startk);
            endk = min(out, endk);
            stepk = floor(stepk/2); %smaller grid
            startk = max(endk - length(startk:stepk:endk)*stepk, startk);
        case 'max'
            startk = max(out, startk);
            endk = min(out + stepk, endk);
            stepk = floor(stepk/2); %smaller grid
            startk = max(endk - length(startk:stepk:endk)*stepk, startk);
    end
    
    if  stepk <=0 || startk > endk
        break;
    end
end

    function [out fval] = bestongrid(stepk, startk, endk)
        %compute error bound for a grid of k
        %best = -Inf;
        out = NaN;
        fval = NaN;
        for i = 1:length(startk:stepk:endk)
            switch lower(type)
                case 'min'
                    k = startk + (i-1)*stepk;
                case 'max'
                    k = endk - (i-1)*stepk;
            end
            fval = feval(k);
            if fval <= 0
                out = k;
                break;
                
            end
        end
    end

    function out = feval(k)
        if ~isnan(fcache(k - startk0 + 1))
            out = fcache(k - startk0 + 1);
        else
            out = fhandle(k);
            fcache(k - startk0 + 1) = out;
        end
    end
end