function out = Rstar(d, n, e, p, fun)
%finds the rate needed to support a given excess distortion at a given blocklength for coin
%flip source
%d - excess distortion
%n - block length (scalar)
%e - excess probability
%fun - which function to use for calculation

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

%starting points - 'persistent' to make optimization faster
persistent x0;
persistent y0;
slack = .1;
tol = 1e-12;
accuracy = 0; %accuracy: 1 - more accurate, 0 - faster
try

    switch lower(fun)
        case 'shannon'
            %achievability via Shannon's upper bound:
            dglb = [tol; 1];
            dgub = [p; Inf];
            options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set');

            if isempty(x0) || isempty (y0)
                x0 = [h(p) - h(d) 1]; %rate
                y0 = [d; 1]; %delta gamma
            end
            out = Shannon();
        case 'normal'
            %normal approximation
            dlb = tol;
            dub = .5-tol;
            options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set');

            if isempty(x0)
                x0 = [h(p) - h(d) 1];
            end
            out = Normal(); 
        case 'spherecoveringc'
            %converse via sphere covering
            if isempty(x0) || isempty(y0) 
                x0 = [h(p) - h(d) 1];
                y0 = [0 1];
            end            
            out = SphereCoveringC();
        case 'generalc'
            taulb = tol;
            tauub = 2;     
            options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set', 'Display', 'off');

            if isempty(x0) || isempty (y0)
                x0 = [h(p) - h(d) 1 - p]; %rate
                y0 = .5; %tau
            end
            out = GeneralC();
        case 'spherecoveringa'
            %achievability via sphere covering 
            if isempty(x0)
                x0 = [h(p) - h(d) 1]; %rate
                %x0(2) = mean(x0);
            end
            options = optimset('TolX',tol, 'MaxFunEvals', 100);
            out = SphereCoveringA();
        case 'permutationa'
            if isempty(x0)
                x0 = [fsolve(@(x)h(p) - h(x) - R,.11) p];
                x0(2) = mean(x0);
            end
            options = optimset('TolX',tol, 'MaxFunEvals', 100);
            out = DPermutationA();
        case 'marton'
            %converse via Marton's error exponent
            qdlb = [p tol];
            qdub = [1/2 1];
            options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set');
            if isempty(x0) && isempty (y0)
                y0 = [1/2 d]; %q0 %delta0
            end
            out = Marton();
        otherwise
            disp('Error: unknown type.')
    end

catch ME
    fprintf('Error: %s\n', ME.message);
    out = NaN;
end


%--------------------------------------------------------------------------
    function out = Shannon()
        %achievability via Shannon's upper bound:
        out = fzero(@(x)optR(x) - e,x0, options);
        %x0(2) = out + slack;
        x0 = out;
        function out = optR(x)
            R = x;
            [dg, Pe] = fmincon(@optdg, y0,[],[],[],[],dglb, dgub, [], options);
            if abs(Pe - e) < 10*tol
                y0 = dg;
            end
            out = Pe;

            function out = optdg(x)
                delta = x(1);
                g = log2(x(2));
                out = Peub(R, n, d, delta, g, p);
            end
        end

    end

%--------------------------------------------------------------------------
    function out = Normal()
        %normal approximation
        if p == 1/2
            out = 1 - h(d) + log2(n)/2/n;
        else
            V = p*(1-p)*(log2((1-p)/p))^2;
            out = h(p) - h(d) + sqrt(V/n)*Qinv(e);
        end
    end

%--------------------------------------------------------------------------
    function out = SphereCoveringC()
        %converse via sphere covering
        if p == 1/2
            out = fzero(@(x)PelbSphereCoveringFair(x,n,d) - e, x0);
            x0(2) = out+slack;
            return;
        end

        %find rstar:
        rstar = fzero(@(r) vol(r*n)- 1 + e, y0);
        %y0(2) = rstar+slack;
        
        rstar = floor(rstar*n);
        while vol(rstar) < 1 - e
            rstar = rstar + 1;
        end
        while vol(rstar) > 1 - e
            rstar = rstar - 1;
        end

        lambda = (1 - e - binocdf(rstar, n, p))/binopdf(rstar +1, n, p);
        out = rate(d);
        fprintf('SphereCoveringC: n = %i, r = %f, lambda = %f, R = %f\n', n, rstar/n, lambda, out);
        
        function out = vol(r)
            out = binocdf(floor(r), n, p);
        end

        function out = rate(d)
            out = (log2( Csum(n, rstar, 'lb', accuracy) + lambda*C(n, rstar+1, 'lb') ) - log2(Csum(n, floor(d*n), 'ub', accuracy )) )/n;
        end

        function out = PelbSphereCoveringFair(R, n, d)
            %lower bound on P_excess for FAIR coin flip source (approx)
            %R - rate
            %n - block length (scalar)
            %D - distortion
            out =1 - Csum(n, floor(n*d), 'ub', accuracy)/2^(n*(1 - R));         
        end
    end
%--------------------------------------------------------------------------

    function out = GeneralC()
        out = fzero(@(x)Pe(x) - e, x0, options);
       
        function out = Pe(R)
            %general converse
            shift = 5;
            [dummy, f] = fmincon(@opttau, y0,[],[],[],[],taulb, tauub, [], options);
            out = 2^(-f)-shift;
            function out = opttau(tau)
                out = 1 - binocdf( ((R + h(d) + log2(1-p))*n + tau*log2(n))/(log2(1-p) - log2(p)),n, p) - 1/n^tau;
                out = - log2(shift + out);
            end
        end
    end

%--------------------------------------------------------------------------
    function out = SphereCoveringA()
        %achievability via sphere covering 
        
        if p == 1/2
        %FAIR coin flip source
            out = fzero(@(x)logPeubSphereCoveringFair(x,n,d) - log(e), x0);
            x0(2) = out +slack;
        else
        %arbitrary p source
            out = fzero(@(x)PeubSphereCovering(x,n,d) - e, x0, options);
            %x0(2) = out +slack;
            x0 = out;
        end
        
        function out = logPeubSphereCoveringFair(R, n, d)
            if n <= 100
                out =  2^(n*R) *log(1 - Csum(n, floor(n*d), 'lb', accuracy) /(2^n));
            else
                out = -2^(n*(R-1))*Csum(n, floor(n*d), 'lb', accuracy);
                
            end
        end
        
        function out = PeubSphereCovering(R, n, d)
            q = max ( 0, (p-d)/(1-2*d));
            if q == 0 
                q = p;
            end
            M = 2^(n*R);
            out = 0;
            for t = 0:n
                out = out + binopdf(t,n,p)*e(n,t);
            end
            %fprintf('SphereCoveringA: n = %i, R = %f, Peub = %f\n', n, R, out);
            
            function out = e(n,t)
                out = 0;
                if 1 
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
    function out = Marton()
        %converse via Marton's error exponent
        Rmax = h(1/2) - h(d) - tol;
        Rmin = h(p) - h(d);
        
        R = fzero(@(x) optR(x) - e, [Rmin Rmax], options); 

        if R == Rmin
           R = fzero(@(x)epsilonRlessRmin(x) - e, [0 Rmin], options);
        end
        out = R;
        %x0 = out;
        function out = optR(x)
            R = x;

            qlb = fzero(@(x)h(x) - h(d) - R,[p 1/2], options); %R(q, d) > R

            if y0(1) < qlb %y0 = q0 d0
                y0(1) = qlb;
            end

            [qd, minusPe] = fmincon(@optqd, y0,[],[],[],[],qdlb, qdub, [], options);
            out = -minusPe;
            if abs(out - e ) < 10*tol
                y0 = qd; 
            end
            function out = optqd(x)
                q = x(1);
                delta = x(2);
                out = -PelbMarton(R, n, d, q, delta, p);
            end

        end
        
        function out = epsilonRlessRmin(R)
            DR = fzero(@(x) h(p) - h(x) - R, [0 p], options);
            out = (DR - d)/(1 - d);
        end
    end

end
