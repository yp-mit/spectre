function out = Dstar(R, n, e, p, fun)
%finds excess distortion for a given rate at a given blocklength for coin
%flip source
%R - rate
%n - block length (scalar)
%e - excess probability
%fun - what function to use for calculation

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%


%starting points - 'persistent' to make optimization faster
persistent x0;
persistent y0;
slack = .1;
try

    switch lower(fun)
        case 'shannon'
            %achievability via Shannon's upper bound:
            tol = 1e-4;
            dglb = [tol; 1];
            dgub = [.5-tol; Inf];
            options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set');

            if isempty(x0) && isempty (y0)
                x0 = [tol 1/2];
                y0 = [.1; 1]; %delta gamma
            end
            out = DShannon();
        case 'normal'
            %normal approximation
            tol = 1e-4;
            dlb = tol;
            dub = .5-tol;
            options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set');

            if isempty(x0)
                x0 = [fsolve(@(x)h(p) - h(x) - R,.11) 1/2];
            end
            out = DNormal(); 
        case 'spherecoveringc'
            %converse via sphere covering
            if isempty(x0) && isempty(y0) 
                x0 = [fsolve(@(x)h(p) - h(x) - R,.11) p];
                y0 = [0 1];
            end            
            out = DSphereCoveringC();
        case 'spherecoveringa'
            %achievability via sphere covering 
            if isempty(x0)
                x0 = [fsolve(@(x)h(p) - h(x) - R,.11) p];
                x0(2) = mean(x0);
            end
            tol = 1e-3;
            options = optimset('TolX',tol, 'MaxFunEvals', 100);
            out = DSphereCoveringA();
        case 'permutationa'
            if isempty(x0)
                x0 = [fsolve(@(x)h(p) - h(x) - R,.11) p];
                x0(2) = mean(x0);
            end
            tol = 1e-3;
            options = optimset('TolX',tol, 'MaxFunEvals', 100);
            out = DPermutationA();
        case 'marton'
            %converse via Marton's error exponent
            tol = 1e-4;
            qdlb = [p tol];
            qdub = [1/2 1];
            options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set');
            if isempty(x0) && isempty (y0)
                %x0 = [tol 1/2];
                y0 = [.1/2 .1]; %q0 %delta0
            end
            out = DMarton();
        otherwise
            disp('Error: unknown type.')
    end

catch ME
    fprintf('Error: %s\n', ME.message);
    out = NaN;
end


%--------------------------------------------------------------------------
    function out = DShannon()
        %achievability via Shannon's upper bound:
        out = fzero(@optD,x0);
        x0(2) = out + slack;
        function out = optD(x)
            D = x;
            [y0, Pe] = fmincon(@optdg, y0,[],[],[],[],dglb, dgub, [], options);
            out = Pe - e;

            function out = optdg(x)
                d = x(1);
                g = log2(x(2));
                out = Peub(R, n, D, d, g, p);
            end
        end

    end

%--------------------------------------------------------------------------
    function out = DNormal()
        %normal approximation
        out = fzero(@(x)Rnormal(x) - R,x0);
        x0 = out + slack;
        function out = Rnormal(D)
            if p == 1/2
                out = 1 - h(D) + log2(n)/2/n;
            else
                V = p*(1-p)*(log2((1-p)/p))^2;
                out = h(p) - h(D) + sqrt(V/n)*Qinv(e);
            end
        end
    end

%--------------------------------------------------------------------------
    function out = DSphereCoveringC()
        %converse via sphere covering
        if p == 1/2
            out = fzero(@(x)PelbSphereCovering(R,n,x) - e, x0);
            x0(2) = out+slack;
            return;
        end

        %find rstar:
        rstar = fzero(@(r) vol(r*n)- 1 + e, y0);
        y0(2) = rstar;
        
        rstar = floor(rstar*n);
        while vol(rstar) > 1 - e
            rstar = rstar - 1;
        end

        lambda = (1 - e - binocdf(rstar, n, p))/binopdf(rstar +1, n, p);
        out = fzero(@(x) R - rate(x), x0);
        x0(2) = out;
        fprintf('SphereCoveringC: n = %i, r = %i, lambda = %f, d = %f\n', n, rstar, lambda, out);
        
        function out = vol(r)
            out = binocdf(floor(r), n, p);
        end

        function out = binosum(l, L, UBorLB)
            out = 0;
            for j = 0:L
                out = out + C(l,j, UBorLB);
            end
        end

        function out = rate(d)
            out = (log2(binosum(n, rstar, 'lb') + lambda*C(n, rstar+1, 'lb')) - log2(binosum(n, floor(d*n), 'ub' )) )/n;
        end

    end
%--------------------------------------------------------------------------
    function out = DSphereCoveringA()
        %achievability via sphere covering 
        
        if p == 1/2
        %FAIR coin flip source
            out = fzero(@(x)logPeubSphereCoveringFair(R,n,x) - log(e), x0);
            x0(2) = out +slack;
        else
        %arbitrary p source
            out = fzero(@(x)PeubSphereCovering(R,n,x) - e, x0, options);
            x0(2) = out +slack;
        end
        
        function out = logPeubSphereCoveringFair(R, n, d)
            if n <= 100
                out =  2^(n*R) *log(1 - binosum(n, floor(n*d)) /(2^n));
            else
                out = -2^(n*(R-1))*binosum(n, floor(n*d));
                
            end
            function out = binosum(l, L)
                out = 0;
                for j = 0:L
                    out = out + nchoosek(l,j);
                end
            end
        end
        
        function out = PeubSphereCovering(R, n, d)
            q = max ( 0, (p-d)/(1-2*d));
            if q == 0 q = p;end
            M = 2^(n*R);
            out = 0;
            for t = 0:n
                out = out + binopdf(t,n,p)*e(n,t);
            end
            fprintf('SphereCoveringA: n = %i, d = %f, Peub = %f\n', n, d, out);
            
            function out = e(n,t)
                out = 0;
                if 0 
                    for T = 0:n
                        out = out + binopdf(T, n, q)*Pwithin(n, t, T);
                    end
                else
                   T = round(q*n); %#ones
                   out = Pwithin(n,t, T);
                end
       
                if n <= 50
                    out = ( 1 - out)^M;
                else
                    out = exp(-M*out);
                end
                
            end
            
            function out = Pwithin(n,t, T)
                %y of type T is within a given x of type t - almost exact lower bound               
                
                out = 0;
                start = max ( 0, ceil( (t+T - n*d)/2) ); %# ones  in 11111 zone (T ones in codewords)
                for k = start:start; %t
                    out = out + C(T, k, 'lb')*C(n-T, t - k, 'lb'); 
                end
  %              fprintf(' d = %f, t = %i, T = %i, #inside d-ball = %f\n', d, t, T, out);
                out = out / C(n, t, 'ub');
            end
        end
    end
 
%--------------------------------------------------------------------------
    function out = DPermutationA()
        %achievability via permutation coding
        
        %arbitrary p source <1/2
        out = fzero(@(x)PeubPermutation(R,n,x) - e, x0, options);
        x0(2) = out;
        
        function out = PeubPermutation(R, n, d)
            q = max ( 0, (p-d)/(1-2*d));
            t = ceil(q*n);
            targetM = 2^(n*R);
            realM = C(n,t,'ub');
            for k = 1:t
                addM = C(n,t-k,'ub') + C(n,t+k,'ub');
                if realM + addM <= targetM
                    realM = realM + addM;
                else
                    break;
                end
            end
            
            out = binocdf(t-k-d*n, n, p) + 1 -  binocdf(t+k+d*n, n, p);
        end
    end

%--------------------------------------------------------------------------
    function out = DMarton()
        %converse via Marton's error exponent
        dmax = fsolve(@(x)h(1/2) - h(x) - R,.11);
        dmin = fsolve(@(x)h(p) - h(x) - R,.11);
        out = fminbnd(@(x) abs(optD(x) - e), dmin, dmax);
        %x0 = out;
        function out = optD(x)
            D = x;

            qlb = fsolve(@(x)h(x) - h(D) - R,.11); %R(q, d) > R

            if y0(1) < qlb
                y0(1) = qlb;
            end

            [y0, Pe] = fmincon(@optqd, y0,[],[],[],[],qdlb, qdub, [], options);
            out = -Pe;

            function out = optqd(x)
                q = x(1);
                delta = x(2);
                out = -PelbMarton(R, n, D, q, delta, p);
            end

        end
    end

end
