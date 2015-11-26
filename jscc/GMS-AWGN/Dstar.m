function out = Dstar(R, n, e, P, fun)
%finds excess distortion for a given rate at a given blocklength for
%transmission of a GMS over an AWGN (noise variance 1)
%R - rate
%n - block length (scalar)
%e - excess probability
%P - channel SNR
%fun - which function to use for calculation

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

%starting points - 'persistent' to make optimization faster
persistent x0;
persistent y0;
persistent nprev;

if isempty(nprev)
    nprev = 0;
end

tol = 1e-12;

options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Display', 'off');
dbar = 1/(1+P);
V1 = 2*dbar^2; %no coding (amplifiers)
V = 2*dbar^2*(2 - dbar^2); %optimal blocklength-n code

Vc = P*(P+2)/2/(P+1)^2;
Vs = .5;

deltan = n - nprev; %distance between consecutive n's
x0correction = sqrt(V1/n)*Qinv(e)/2/n^(3/2)*deltan;
y0slack = .1;

k = ceil(R*n);

if ~isempty(x0)
    x0 = x0 - x0correction;
    d0 = x0;
else
    d0 = 1.8*dbar;
end

switch lower(fun)
    case 'aseparate'
        out = ASeparate();
    case 'ajoint'
        %Achievability via JSCC coding
        options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set', 'Display', 'off');
        out = AJoint();
    case 'cht'
        %Converse via hypothesis testing
        out = Cht();
    case 'nocoding'
        %no coding
        out = ANoCoding();
    case 'approxnocoding'
        out = ApproxNoCoding();
    case 'approx'
        %gaussian approximation
        out = Approx();
    otherwise
        disp('Error: unknown type.')
end



%--------------------------------------------------------------------------
    function out = ANoCoding()
        %No coding
        
        out = fzero(@(d)1 - e - chi2cdf(n*d/dbar, n), d0);
    end

%--------------------------------------------------------------------------
    function out = ApproxNoCoding()
        %No coding
        out = dbar + sqrt(V1/n)*Qinv(e);
    end
%--------------------------------------------------------------------------
    function out = Approx()
        %No coding
        %C = .5*log(1 + P);
       % Rd = .5*log(1/d);
        %out = fzero(@(d)C + .5*log(d) - sqrt(1/n*(Vc + Vs))*Qinv(e), d0);
        out = dbar + sqrt(V/n)*Qinv(e);
    end


%--------------------------------------------------------------------------
    function out = Cht()
        %Converse via meta-converse
        
        tau = thres(k);
        
        out = fzero(@(d)Qrob(d, tau) - 1, d0);
        updatex0(out);
              
        
        function out = Qrob(d, tau)
            
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
            %out = fzero(@(tau)Prob(tau) - 1 + e, tau0, options);
            out = fzero(@(tau)Prob(tau) - 1 + e, tau0, options);
            
            
            
            function out = Prob(tau)
                
                if n < 200
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
        %dP/dQ is bounded by a constant. Find that constant.
        [~, gamma2] = fminbnd(@(x)-ncx2pdf(n*x, n, n*P)./chi2pdf(n*x/(1+P), n)*(1+P), 0, 10);
        gamma2 = -gamma2;

        out = fzero(@(d)bound(d) - e, d0);
        updatex0(out);
%        fprintf('AJoint: n = %i, R = %f \n', n, out);
        
        
        function out = bound(d)           
            %optimize with respect to r and g; y0 = [r gg] - TOO SLOW
            if d > 1
                out = NaN;
                return;
            end
            out = goal([d 3]);
               
            
            function out = goal(dg)
                d = dg(1);
                    
                gamma1 = log(n)/2*dg(2);
                
                r0 = sqrt(1 - d);
                a = max(0, r0 - sqrt(d));
                b = r0 + sqrt(d);
                
                
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
                    
                    %handle small distances from the origin
                    out(x <= sqrt(d) - r0) = -Inf;
                end
            end
        end
    end


%--------------------------------------------------------------------------
    function out = ASeparate()
        %achievability via separate coding
        %R = 1
        
        if isempty(y0)
            %first pass
            startL = 1;
            endL = n;
        else
            startL = max( ceil((y0 - 4*y0slack)*n), 1);
            endL = min( ceil((y0 + 4*y0slack)*n), n );
        end
        
        %precompute channel error probability
        PeCh = ones(1, endL);
        for L = startL:endL %exp(r) - the number of messages passed from the source enc to the channel enc
            PeCh(L) = PeubChannel(P, n, L);
        end
        
        Peubseparate(k, d0, PeCh, e);
        out = fzero(@(d)Peubseparate(k, d, PeCh, e) - e, d0);
        updatex0(out);
        [err, curL] = Peubseparate(k, out, PeCh, e);
        updatey0(curL/n); 
    end

%--------------------------------------------------------------------------
%Utilities for faster computation
    function updatex0(x)
        if ~isnan(x) && n > 100
            x0 = x;
            nprev = n;
        end
    end

%--------------------------------------------------------------------------
%Utilities for faster computation
    function updatey0(y)
        if ~isnan(y) && n > 100
            y0 = y;
        end
    end
end


