function out = Dstar(R, n, e, fun)
%finds excess distortion for a given rate at a given blocklength for
%Gaussian source with unit variance
%R - rate
%n - block length (scalar)
%e - excess probability
%fun - which function to use for calculation

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%


%starting points - 'persistent' to make optimization faster
persistent x0;
persistent y0;

try

    switch lower(fun)
        case 'shannon'
            %achievability via Shannon's upper bound:
            tol = 1e-4;
            qglb = [2^(-2*R); 1];
            qgub = [1; Inf];
            options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set');

            if isempty(x0) && isempty (y0)
                x0 = [2^(-2*R) 1];
                y0 = [(2^(-2*R)+1)/2 1]; %q gamma
            end
            out = DShannon();
        case 'normal'
            %normal approximation
            out = DNormal(); 
        case 'spherecoveringc'
            %converse via sphere covering
            out = DSphereCovering();
        case 'spherecoveringa'
            %achievability via sphere covering
            if isempty(x0) && isempty (y0)
                x0 = [2^(-2*R) .5]; %distortion
             %   y0 = [0 1];        %auxiliary parameter q
            end 
            out = DSphereCoveringA();
        otherwise
            disp('Unknown type.')
    end

catch ME
    fprintf('Error: %s\n', ME.message);
    out = NaN;
end


%--------------------------------------------------------------------------
    function out = DShannon()
        %achievability via Shannon's upper bound:
        out = fzero(@(x)Pexcess(x) - e,x0);
        x0 = out;
        function out = Pexcess(x)
            D = x;
            [y0, Pe] = fmincon(@optqg, y0,[],[],[],[],qglb, qgub, [], options);
            out = Pe;

            function out = optqg(y)
                q = y(1);
                g = log2(y(2));
                out = PeubShannon(R, n, D, q, g);
            end
        end

    end

%--------------------------------------------------------------------------
    function out = DNormal()
        out = 2^(-2*R)*(1 + sqrt(2/n)*Qinv(e));
    end

%--------------------------------------------------------------------------
    function out = DSphereCovering()
        %converse via sphere covering

        %find rstar:
        factor = 10;
        rstar = fzero(@(x)chi2cdf(x^2*n, n) - chi2cdf(factor*n, n) + e, 1);
        out = fminbnd(@(x) abs(R - rate(x)),2^(-2*R),1);


        function out = rate(x)
            out = 1/2*log2(rstar^2/x);
        end
    end

%--------------------------------------------------------------------------
    function out = DSphereCoveringA()
        %achievability via sphere covering
        M = 2.^(n*R);
        out = fzero(@(x)Pexcess(x) - e, x0);
    function out = Pexcess(x)
        d = x;
        s0 = 1; %variance of source
        r0 = sqrt(s0^2 - d); %distance from the center of the representation sphere
        
        out = quad(@density, r0 - sqrt(d), r0 + sqrt(d))...
            + 1 - chi2cdf((r0 + sqrt(d))^2*n,n)...
            + chi2cdf((r0 - sqrt(d))^2*n,n);
        out = real(out);
        
        function out = density(r)
          cosa = (r0^2 + r.^2 - d)./(2*r0.*r);
          sina = (1 - cosa.^2).^(1/2);
          out = ub(sina).*chi2pdf(n*r.^2, n).*n*2.*r;         
        end
        
        function out = ub(sina)
            if (n < 50)
                A = 1/sqrt(pi).*gamma(n/2+1)./n./gamma((n-1)/2+1);
                out = (1 - A.*(sina).^(n-1)).^M;
            else
                out = exp(-M*sina.^(n-1)./sqrt(2*pi*n));
            end
        end
    end  
        
    end



end
