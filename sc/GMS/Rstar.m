function out = Rstar(d, n, e, fun)
%finds rate that achieves a given excess distortion at a given blocklength for
%Gaussian source with unit variance
%d - excess distortion
%n - block length (scalar)
%e - excess probability
%fun - which function to use for calculation

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

%starting points - 'persistent' to make optimization faster
persistent x0;
%rate distortion
Rd = -1/2*log2(d);
%precision options
tol = 1e-10;
options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set', 'Display', 'off');

switch lower(fun)
    case {'shannon', 'spherecoveringa', 'generalc', 'generalcopt'}
        if isempty(x0)
            x0 = [Rd 3*Rd]; %rate
        end
        out = Generic();
    case 'normal'
        %normal approximation:
        out = Normal();
    case 'spherecoveringc'
        %converse via sphere covering
        out = SphereCoveringC();
    case 'rogersa'
        out = RogersA();
    otherwise
        disp('Unknown type.')
end


%--------------------------------------------------------------------------
    function out = Generic()
        out = fzero(@(x)goal(x) - e, x0, options);
        function out = goal(x)
            out = Pexcess(x, d, n, fun);
        end
    end

%--------------------------------------------------------------------------
    function out = Normal()
        %Normal approximation
        out = - .5*log2(d) + 1/sqrt(2*n)*Qinv(e)*log2(exp(1));
    end

%--------------------------------------------------------------------------
    function out = SphereCoveringC()
        %converse via sphere covering
        
        %find rstar:
        factor = 3;
        rstar = fzero(@(x)chi2cdf(x^2*n, n) - chi2cdf(factor*n, n) + e, 1);
        out = rate(d);
        
        function out = rate(x)
            out = 1/2*log2(rstar^2/x);
        end
    end

%--------------------------------------------------------------------------
    function out = RogersA()
        rstar = fzero(@(x)1 - chi2cdf(x^2*n, n) - e, 1);
        out = rogers(rstar/sqrt(d));
        function out = rogers(r)
            if r >= n
                out = log2(exp(1))+ log2((n*log(n)+n*log(log(n)) + 5*n)) + n*log2(r);
            elseif r >= n/log(n)
                out = log2(n*(n*log(n)+n*log(log(n)) + 5*n)) + n*log2(r);
            elseif r >2
                out = log2(7^(4*log(7)/7)/4*sqrt(2*pi)) + log2(n^(3/2)*((n-1)*log(r*n)+(n-1)*log(log(n)+log(n)/2+log((pi*sqrt(2*n))/(sqrt(pi*n)-2))))) - log2(r*(1-2/log(n))*(1 - 2/sqrt(pi*n))*(log(n))^2) + n*log2(r);
            else
                out = log2(sqrt(2*pi)) + log2(sqrt(n)*((n-1)*log(r*n)+(n-1)*log(log(n)+log(n)/2+log((pi*sqrt(2*n))/(sqrt(pi*n)-2))))) - log2(r*(1-2/log(n))*(1 - 2/sqrt(pi*n))) + n*log2(r);
            end
            out = out/n;
        end
    end

end
