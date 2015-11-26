function out = Pexcess(R, d, n, fun)
%finds excess probability at a given rate, distortion and blocklength for
%Gaussian source with unit variance
%R - rate
%d - excess distortion
%n - block length (scalar)
%outputs excess probability
%fun - which function to use for calculation

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

%starting points - 'persistent' to make optimization faster

persistent y0;
%rate distortion
Rd = -1/2*log2(d);
%precision options
tol = 1e-10;
options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set', 'Display', 'off');

switch lower(fun)
    case 'shannon'
        %achievability via Shannon's upper bound:
        qglb = [tol; 1];
        qgub = [1; Inf];
        if isempty (y0)
            y0 = [d 1]; %q gamma
        end
        out = Shannon();
    case 'normal'
        %normal approximation:
        out = Normal();
    case 'spherecoveringc'
        %converse via sphere covering
        out = SphereCoveringC();
    case 'generalc'
        taulb = tol;
        tauub = 2;
        if isempty (y0)
            y0 = .5; %tau
        end
        out = GeneralC();
    case 'generalcopt'
        y0lb = [tol tol];
        y0ub = [1 2];
        if isempty (y0)
            y0 = [d .5]; % dd tau
        end
        out = GeneralCOpt();
    case 'spherecoveringa'
        %achievability via sphere covering
        out = SphereCoveringA();
    case 'rogersa'
        out = RogersA();
    otherwise
        disp('Unknown type.')
end


%--------------------------------------------------------------------------
    function out = Shannon()
        %achievability via Shannon's upper bound:
        [qg, Pe] = fmincon(@optqg, y0,[],[],[],[],qglb, qgub, [], options);
        if Pe < .5 && Pe > 0
            y0 = qg;
        end
        out = Pe;
        
        function out = optqg(y)
            q = y(1);
            g = log2(y(2));
            out = PeubShannon(R, n, d, q, g);
        end
    end

%--------------------------------------------------------------------------
    function out = Normal()
        %normal approximation
        out = Q((R + log2(d)/2)*sqrt(2*n)/log2(exp(1)));
    end

%--------------------------------------------------------------------------
    function out = SphereCoveringC()
        %converse via sphere covering
        
        %find rstar:
        r = sqrt(d)*2^R;
        out = 1 - chi2cdf(r^2*n, n);
    end

%--------------------------------------------------------------------------
    function out = GeneralC()
        %general converse
        shift = 5;
        [y0, f] = fmincon(@opttau, y0,[],[],[],[],taulb, tauub, [], options);
        out = 2^(-f)-shift;
        function out = opttau(tau)
            out = 1 - chi2cdf( (1 - log(1/d) + 2*R/log2(exp(1)))*n + 2*tau*log(n),n) - 1/n^tau;
            out = - log2(shift + out);
        end
    end

%--------------------------------------------------------------------------
    function out = GeneralCOpt()
        %general converse optimized over dd
        shift = 5;
        [y0, f] = fmincon(@goal, y0,[],[],[],[],y0lb, y0ub, [], options);
        out = 2^(-f)-shift;
        function out = goal(y)
            dd = y(1);
            tau = y(2);
            out = 1 - chi2cdf( (d/dd - log(1/dd) + 2*R/log2(exp(1)))*n + 2*tau*log(n),n) - 1/n^tau;
            out = - log2(shift + out);
        end
    end

%--------------------------------------------------------------------------
    function out = SphereCoveringA()
        %achievability via sphere covering
        s0 = 1; %variance of source
        r0 = sqrt(s0^2 - d); %distance from the center of the representation sphere
        
        out = quad(@density, r0 - sqrt(d), r0 + sqrt(d), tol)...
            + 1 - chi2cdf((r0 + sqrt(d))^2*n,n)...
            + chi2cdf((r0 - sqrt(d))^2*n,n);
        out = real(out);
        
        function out = density(r)
            cosa = (r0^2 + r.^2 - d)./(2*r0.*r);
            sina = (1 - cosa.^2).^(1/2);
            out = ub(sina).*chi2pdf(n*r.^2, n).*n*2.*r;
        end
        
        function out = ub(sina)
            if (n < 20)
                M = 2.^(n.*R);
                A = 1/sqrt(pi).*gamma(n/2+1)./n./gamma((n-1)/2+1);
                out = (1 - A.*(sina).^(n-1)).^M;
            else
                out = exp(-2.^(n.*R+(n-1).*log2(sina))./sqrt(2*pi*n));
            end
        end
    end

%--------------------------------------------------------------------------
    function out = RogersA()
        r = fzero(@(x)R - rogers(abs(x)/sqrt(d)), sqrt(d));
        out = 1 - chi2cdf(r^2*n, n);
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
