function out = Rstar(P, d, e, n, fun)
%finds the rate required for a given excess distortion at a given blocklength for
%transmission of a GMS over an AWGN (noise variance 1)
%d - excess distortion (d/sigma^2)
%P - allowable power
%e - excess probability
%n - block length (scalar)
%fun - which function to use for calculation

%EVERYTHING IN NATS

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

%starting points - 'persistent' to make optimization faster
persistent x0;
persistent y0;
persistent x0prev;
persistent y0prev;
lenprev = 10;

tol = 1e-12;

Capacity = .5*log(1 + P);
Rd = 1/2*log(1/d);
Rlim = Capacity/Rd;
Vc = P*(P+2)/2/(P+1)^2;
Vs = .5;
deltan = 10; %distance between consecutive n's
x0correction = 1/Rd*sqrt((Vc + Rlim*Vs))*Qinv(e)/2/n^(3/2)*deltan;
x0slack = 2*x0correction;
y0slack = .1*Capacity;

%if length(x0prev) > 1
%    x0slack = min(.5, 6*std(x0prev));
%end
if length(y0prev) > 1
    y0slack = min(.5, 6*std(y0prev));
end



options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set', 'Display', 'off');

switch lower(fun)
    case 'approx'
        %gaussian approximation
        out = Approx();
    case 'cht'
        %Converse via list decoding (hypothesis testing)
        out = Cht();
    case 'cgeneral'
        taulb = tol;
        tauub = 10;
        if isempty(x0)
            x0 = Rlim; %upper bound for rate
            y0 = 1; %tau
        end
        
        out = Cgeneral();
    case 'approxseparate'
        %gaussian approximation
        out = ApproxSeparate();
    case 'aseparate'
        %Achievability via separate coding
        out = ASeparate();
    case 'ajoint'
        %Achivability via JSCC coding
        if isempty(x0) || isempty(y0)
            x0 = Rlim; %upper bound for rate
            y0 = 20; %gg
        end
        out = AJoint();
    otherwise
        disp('Error: unknown type.')
end




%--------------------------------------------------------------------------
    function out = Approx()
        %Gaussian approximation
        out = Rlim - 1/Rd*sqrt((Vc + Rlim*Vs)/n)*Qinv(e);
    end

%--------------------------------------------------------------------------
    function out = ApproxSeparate()
        %Gaussian approximation of separate source-channel coding
        [~, opt] = fmincon(@goal,[e/2, e/2],[1 1], e,[],[], [0 0], [e e], [], options);
        out = Rlim - 1/sqrt(n)/Rd*opt;
        
        function out = goal(arg)
            out = sqrt(Vs*Rlim)*Qinv(arg(1)) + sqrt(Vc)*Qinv(arg(2));
        end
    end
%--------------------------------------------------------------------------

    function out = Cgeneral()
        if x0 >= Rlim
            %first pass
            startk = 1;
            endk = floor(Rlim*n);
        else
            startk = floor( max( 1, (x0-x0slack)*n ) );
            endk = floor( min( Rlim*n, (x0+x0slack)*n ) );
        end
        
        out = fzerointeger(@(x)Pe(x/n) - e, startk, endk, 'max');
        out = out/n;
        updatex0(out);
        
        %out = fzero(@(x)Pe(x) - e, x0, options);
        
        function out = Pe(R)
            shift = 10;
            [~, f] = fmincon(@(x)-logPe(x), y0,[],[],[],[],taulb, tauub, [], options);
            %           tau = 2;
            %           out = logPe(tau);
            out = exp(-f) - shift;
            function out = logPe(tau)
                gamma = log(n)/2*tau;
                k = round(R*n);
                
                if n < 500
                    b = 10*k/n; %Ev0 = k/n
                else
                    b = 2*k/n;
                end
                
                out = quad(@density, 0, b) - exp(-gamma);
                out = real(log(out + shift));
                function out = density(v0)
                    % v - noncentral chi square (n) /n
                    % v0 - central chi square (k) /n
                    
                    out = (1 - ncx2cdf(thres(v0),n,n/P)).*chi2pdf(n*v0,k)*n;
                    
                    function out = thres(v0)
                        out = n/2*(log(1+P)+1) + k/2*(1 - log(1/d));
                        out = 2*(1+P)/P*(out + gamma) - (1+P)/P*n*v0;
                    end
                end
            end
        end
    end

%--------------------------------------------------------------------------

    function out = Cht()
        if isempty(x0)
            %first pass
            startk = 1;
            endk = floor(Rlim*n);
        else
            x0 = x0 + x0correction;
            startk = floor( max( 1, (x0-x0slack)*n ) );
            endk = floor( min( Rlim*n, (x0+x0slack)*n ) );
        end
        
        out = fzerointeger(@(k)Qrob(k) - 1, startk, endk, 'max');
        out = out/n;
        updatex0(out);
        
        function out = Qrob(k)
            tau = thres(k);
            b = sqrt( tau*n/k/d );
            out = quad(@density, 0, b);
            function out = density(r)
                % v - noncentral chi square (n) /n
                % v0 - central chi square (k) /n
                
                out = ncx2cdf(n*(tau - k/n*d*r.^2)/P, n, n*(1 + 1/P)).*r.^(k-1).*k;
                
            end
        end
        
        function out = thres(k)
            %find tau:
            tau0 = 1 + k/n + 2*sqrt( (n*Vc + k*Vs)/n^2 )*Qinv(e);
            out = fzero(@(tau)Prob(tau) - 1 + e, tau0, options);
            
            
            
            function out = Prob(tau)
                
                if n < 500
                    b = 10*k/n; %Ev0 = k/n
                else
                    b = 2*k/n;
                end
                
                out = quad(@density, 0, b);
                function out = density(v0)
                    % v - noncentral chi square (n) /n
                    % v0 - central chi square (k) /n
                    
                    out = ncx2cdf(n*(tau - v0)*(1+P)/P, n, n/P).*chi2pdf(n*v0,k)*n;
                    
                end
            end
        end
    end
%--------------------------------------------------------------------------
    function out = AJoint()
        if x0 >= Rlim
            %first pass
            startk = 1;
            endk = floor(Rlim*n);
        else
            startk = max(1, floor((x0-x0slack)*n ));
            endk = min(floor(Rlim*n), floor((x0+x0slack)*n));
        end
        
        r0 = sqrt(1 - d);
        a = r0 - sqrt(d);
        b = r0 + sqrt(d);
        
        %dP/dQ is bounded by a constant. Find that constant.
        [~, gamma2] = fminbnd(@(x)-ncx2pdf(n*x, n, n*P)./chi2pdf(n*x/(1+P), n)*(1+P), 0, 10);
        gamma2 = -gamma2;
        
        boundCache = NaN(1, floor(Rlim*n));
        out = fzerointeger(@(k)bound(k) - e, startk, endk, 'max');
        out = out/n;
        updatex0(out);
        %        fprintf('AJoint: n = %i, R = %f \n', n, out);
        
        
        function out = bound(k)
            if ~isnan(boundCache(k))
                out = boundCache(k);
                return;
            end
            
            %optimize with respect to r and g; y0 = [r gg] - TOO SLOW
            [g, out] = fmincon(@goal,y0,[],[],[],[], 0, 20, [], options);
            if ~isnan(g)
                y0 = g;
            end
            
            boundCache(k) = out;
            
            function out = goal(g)
                gamma1 = log(n)/2*g;
                
                A = n/2*log(1+P) + n/2 - log(gamma1) - log(gamma2);
                
                
                Fouter = @(v0)chi2pdf(k*v0, k).*k.*Iinner(v0);
                Ftail = @(v0)(1 - ncx2cdf(thresv(v0)*n,n,n/P)).*chi2pdf(k*v0,k).*k;
                %   Ftail = @(v0) exp(n/P*t./(1-2*t))./(1-2*t).^(n/2).*exp(-t*thresv(v0)).*chi2pdf(k*v0,k).*k;
                out = quad(Fouter, a^2 + tol, b^2- tol)...                      %exponent < 1
                    + quad(Ftail, a^2 + tol, b^2 - tol)...                      %exponent = 1
                    + chi2cdf((a^2 + tol)*k, k) + 1 - chi2cdf((b^2 - tol)*k, k)...      %nontypical source realization
                    + exp(1-gamma1);
                
                %out = quad(Ftail, a^2 + tol, b^2 - tol); %most significant
                %term
                
                function out = Iinner(v0)
                    %exponent < 1
                    out = zeros(1, numel(v0));
                    for j0 = 1:numel(v0)
                        Finner = @(v) ncx2pdf(n*v, n, n/P).*n.*exp(P/2/(1+P)*n*(v - thresv(v0(j0))));
                        out(j0) = quad(Finner, 0, thresv(v0(j0)) );
                    end
                    if any(isnan(out))
                        return;
                    end
                end
                
                function out = thresv(x)
                    %threshold for v as a function of v0
                    out = max( 0, 2*(1+P)/P/n*(A + logPdball(x)) );
                end
                
                
                function out = logPdball(x)
                    %x - distance squared from the origin (x*k central chi square k)
                    %lower bound
                    cosa = (r0^2 + x - d)./(2*r0.*sqrt(x));
                    sina = (1 - cosa.^2).^(1/2);
                    out = loggamma(k/2+1, 'LB') - loggamma((k-1)/2 +1, 'UB')- .5*log(pi)-log(k) + (k-1)*log(sina);
                end
                
            end
        end
    end


%--------------------------------------------------------------------------
    function out = ASeparate()
        %achievability via separate coding
        
        if isempty(x0) || isempty(y0)
            %first pass
            startk = 1;
            endk = floor(Rlim*n);
            
            startL = 1;
            endL = n;
        else
            startk = max( ceil(n*(y0 - y0slack)), 1);
            endk = min( ceil(n*(y0 + y0slack)), floor(Rlim*n) );
            
            startL = max( ceil((x0 - 4*x0slack)*n), 1);
            endL = min( ceil((x0 + 4*x0slack)*n), n );
        end
        
        %precompute channel error probability
        PeCh = ones(1, endL);
        for L = startL:endL %exp(r) - the number of messages passed from the source enc to the channel enc
            PeCh(L) = PeubChannel(P, n, L);
        end
        
        kstar = fzerointeger(@(k)Peubseparate(k, d, PeCh, e) - e, startk, endk, 'max');
        [err, curL] = Peubseparate(kstar, d, PeCh, e);
        out = kstar/n;
        updatey0(out);
        updatex0(curL/n);
    end


%--------------------------------------------------------------------------
%Utilities for faster computation
    function updatex0(x)
        if ~isnan(x)
            x0 = x;
            if length(x0prev) < lenprev
                x0prev = [x0prev x0];
            else
                x0prev(1:lenprev - 1) = x0prev(2:lenprev);
                x0prev(lenprev) = x0;
            end
        end
    end


    function updatey0(y)
        if ~isnan(y)
            y0 = y;
            if length(y0prev) < lenprev
                y0prev = [y0prev y0];
            else
                y0prev(1:lenprev - 1) = y0prev(2:lenprev);
                y0prev(lenprev) = y0;
            end
        end
    end
end
