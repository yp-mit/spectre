function out = Dstar(R, n, e, p, delta, fun)
%finds excess distortion for a given rate at a given blocklength for coin
%flip source transmitted over a BSC
%R - rate
%n - block length (scalar)
%e - excess probability
%p - source bias
%delta - channel crossover probability
%fun - which function to use for calculation

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

%starting points - 'persistent' to make optimization faster
persistent x0;
persistent y0;
persistent nprev;

if isempty(nprev)
    nprev = n;
end

tol = 1e-12;
accuracy = 1; %accuracy: 1 - more accurate, 0 - faster



options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Display', 'off');
dbar = fsolve(@(x)R*(h(p) - h(x)) - 1 + h(delta),.11, options);

Rd = h(p) - h(dbar);

Vc = delta*(1-delta)*(log2((1-delta)/delta))^2;
Vs = p*(1-p)*(log2((1-p)/p))^2;
lambda = log2((1-dbar)/dbar);

deltan = n - nprev; %distance between consecutive n's
x0correction = sqrt(Vc + Vs)/lambda*Qinv(e)/2/n^(3/2)*deltan;
slack = 4*abs(x0correction);
y0slack = .1;
k = ceil(R*n);

if ~isempty(x0)
    x0 = x0 - x0correction;
    d0 = x0;
else
    d0 = 2*dbar;
end

switch lower(fun)
    case 'random'
        %achievability via random coding:
        out = ARandom();
    case 'aseparate'
        out = ASeparate();
    case 'ajoint'
        %Achievability via JSCC coding
        options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set', 'Display', 'off');
        out = AJoint();
    case 'clistdecoding'
        %Converse via hypothesis testing
        out = Cht();
    case 'nocoding'
        %no coding
        out = ANoCoding();
    case 'approx'
        %gaussian approximation
        out = Approx();
    otherwise
        disp('Error: unknown type.')
end



%--------------------------------------------------------------------------
    function out = ANoCoding()
        %No coding fair or biased source
        out = fzero(@(d)1 - e - binocdf(n*d, n, delta), 2*delta);
    end


% %--------------------------------------------------------------------------
    function out = ASeparateFair()
        %achievability via separate coding (fair source)
        
        if isempty(x0) || isempty(y0)
            %first pass
            startr = 0; %r = floor(nd)
            endr = k;
            
            startL = 1;
            endL = n;
        else
            startr = max( ceil(k*(x0 - slack)), 0);
            endr = min( ceil(k*(x0 + slack)), k );
            
            startL = max( ceil((y0 - slack)*n), 1);
            endL = min( ceil((y0 + slack)*n), n );
        end
        
        %precompute channel error probability
        PeCh = ones(1, endL);
        for L = startL:endL %2^L - the number of messages passed from the source enc to the channel enc
            PeCh(L) = PeubChannel(n, delta, L, accuracy);
        end
        
        %compute error bound for all d
        best = 0;
        out = NaN;
        for r = startr:endr
            [cur, curL] = Peubseparate(k, r/k, startL, endL, PeCh, accuracy);
            if (cur <= e) && (cur > best)
                best = cur;
                bestL = curL;
                out = r;
            end
        end
        
        out = out/n;
        
        if ~isnan(out)
            x0 = out;
            y0 = bestL/n;
        end
        
    end

%--------------------------------------------------------------------------
    function out = AJoint()
        if isempty(x0)
            %first pass
            startd = ceil(n*dbar);
            endd = floor(3*n*delta);
        else
            startd = ceil( n*max(x0-slack, dbar) );
            endd = floor( n*min(x0+slack, 3*delta) );
        end
        
        %precompute Binomial coeff
        sizeC = n+1;
        logCub = NaN(sizeC, sizeC);
        for kk = 0:sizeC
            for j = 0:sizeC
                logCub(kk+1, j+1) = logC(kk, j, 'ub');
            end
        end
        
        logClb = NaN(sizeC, sizeC);
        for kk = 0:sizeC
            for j = 0:sizeC
                logClb(kk+1, j+1) = logC(kk, j, 'lb');
            end
        end
        
        
        boundCache = NaN(1, n );
        out = fzerointeger(@(d)bound(d) - e, startd, endd, 'min');
        out = out/n;
        updatex0(out);
        %fprintf('AJoint: n = %i, d = %f \n', n, out);
        
        function out = bound(d)
            if ~isnan(boundCache(d))
                out = boundCache(d);
                return;
            end
            
            
            %optimize with respect to r and g; y0 = [r gg]
            %out = goal(y0);
            
            g0 = 10; %gg
            
            [gg, out] = fmincon(@goal,g0,[],[],[],[], 0, Inf, [], options);
            
            boundCache(d) = out;
            
            function out = goal(gg)
                
                %probability of error
                g = .5*log(n)*gg; %natural log
                
                out = 0;
                if p == 1/2
                    rholb = Pdballfair(k, d);
                    for t = 0:n
                        out = out + binopdfub(t, n, delta)...
                            *2^(- max(0, n*log2(2 - 2*delta) - t*log2(1-delta) + t*log2(delta) + log2(rholb) - log2(g)) );
                    end
                else
                    for s = 0:k
                        rholb = rho(k, s);
                        
                        err = 0;
                        for t = 0:n
                            err = err + binopdfub(t, n, delta)...
                                *2^(- max(0, n*log2(2 - 2*delta) - t*log2(1-delta) + t*log2(delta) + log2(rholb) - log2(g)) );
                        end
                        out = out + binopdfub(s, k, p)*err;
                    end
                end
                out = out + exp(1-g);
                
                
                function out = logCcache(n, k, UBorLB)
                    if n < k || k < 0
                        out = 0;
                    else
                        switch lower(UBorLB)
                            case 'ub'
                                out = logCub(n+1, k+1);
                            case 'lb'
                                out = logClb(n+1, k+1);
                        end
                    end
                end
                
                function out = rho(k,s)
                    d1 = d/n;
                    T = ceil((p-d1)/(1-2*d1)*k);
                    j0 = ceil( (T + s - k*d1)/2 );
                    out = 0;
                    for j = j0:s
                        out = out + 2^(logCcache(T, j, 'lb') + logCcache(k-T, s-j, 'lb') - logCcache(k, s, 'ub'));
                    end
                end
                
                function out = Pdballfair(k, d)
                    out = 2^(-k)*Csum(k, d, 'lb', accuracy);
                end
                
                function out = binopdfub(k, n, p)
                    out = 2^(logCcache(n, k, 'ub') + k*log2(p) + (n - k)*log2(1 - p));
                end
            end
        end
    end



%--------------------------------------------------------------------------
    function out = ASeparate()
        %achievability via separate coding
        if isempty(y0) || n <= 100
            startL = 1;
            endL = n;
        else
            startL = max( ceil((y0 - y0slack)*n), 1);
            endL = min( ceil((y0 + y0slack)*n), n );
        end
        
        d0  = fzero(@(d)Peubseparate(k, n, p, delta, d, startL, endL, -1, options) - e, d0, options);
        
        [out, Lstar] = Peubseparate(k, n, p, delta, d0, startL, endL, -1, options);
        updatey0(Lstar/n);
        
        if 1
            %attempt to improve upon the normal approximation solution
            startd = ceil( k*max(d0 - .05, dbar) );
            endd = floor( k*min(d0, p) );
            
            Pech = PeubChannel(n, delta, Lstar, accuracy);
            out = fzerointeger(@(dd)PeubSource(k, p, dd/k, Lstar, accuracy) + Pech - e, startd, endd, 'min');
            out = out/k;
        else
            %verify that normal approximation solution (d, Lstar) leads to an
            %error < e
            if PeubSource(k, p, d0, Lstar, accuracy) + PeubChannel(n, delta, Lstar, accuracy) <= e
                out = d0;
            else
                out = NaN;
            end
        end
        updatex0(out);
    end

%--------------------------------------------------------------------------
    function out = Cht()
        %Converse via meta-converse
        if isempty(d0)
            %first pass
            startd = ceil(n*dbar);
            endd = floor(3*n*delta);
        else
            startd = ceil( n*max(d0-slack, dbar) );
            endd = floor( n*min(d0+slack, 3*delta) );
        end
        
        
        %precompute Binomial coeff
        sizeC = n+1;
        logCub = NaN(sizeC, sizeC);
        for k = 0:sizeC
            for j = 0:sizeC
                logCub(k+1, j+1) = logC(k, j, 'ub');
            end
        end
        
        logClb = NaN(sizeC, sizeC);
        for k = 0:sizeC
            for j = 0:sizeC
                logClb(k+1, j+1) = logC(k, j, 'lb');
            end
        end
        
        %possible values of r that result in different probabilities
        r = NaN(n, n);
        lp = log2((1-p)/p);
        ldelta = log2((1 - delta)/delta);
        for k = 0:n
            for j = 0:n
                r(k+1, j+1) = k*lp + j*ldelta;
            end
        end
        
        r = unique(sort(r(:)));
        
        boundCache = NaN(1, n);
        b = beta(n, e);
        out = fzerointeger(@(d)bound(d, b), startd, endd, 'min');
        out = out/n;
        updatex0(out);
        %fprintf('CListDecoding: n = %i, d = %f \n', n, out);
        
        
        
        function out = beta(n, e)
            k = n;
            %find rstar and lambda:
            Pc = 0;
            prev = -Inf;
            %Compute rstar and lambda
            
            for i = 1:length(r)
                if  r(i) - prev < tol
                    %to exclude dublicates - unique doesn't work well with
                    %reals
                    continue;
                end
                
                add = 0;
                for T = 0:k
                    for t = 0:n
                        if  T*lp + t*ldelta > prev && T*lp + t*ldelta <= r(i)
                            %                             if T*lp + t*ldelta >  k*p*lp + n*delta*ldelta
                            %                              T*lp + t*ldelta
                            %                              2^( log2(Ccache(n, t, 'ub')) + t*log2(delta) + (n-t)*log2(1-delta) + log2(Ccache(k, T, 'ub')) + T*log2(p) + (k-T)*log2(1-p))
                            %                             end
                            add = add + 2^( logCcache(n, t, 'ub') + t*log2(delta) + (n-t)*log2(1-delta) + logCcache(k, T, 'ub') + T*log2(p) + (k-T)*log2(1-p));
                        end
                    end
                end
                
                if Pc + add < 1 - e
                    Pc = Pc + add;
                else
                    lambda = (1 - e - Pc)/add;
                    rstar = prev;
                    break;
                end
                prev = r(i);
            end
            
            %Compute beta_(1-e)
            out = 0;
            for T = 0:k
                for t = 0:n
                    if  T*lp + t*ldelta < rstar
                        out = out + 2^(logCcache(n, t, 'lb')+logCcache(k, T, 'lb') - 2*n);
                    elseif T*lp + t*ldelta - rstar < tol
                        out = out + lambda*2^(logCcache(n, t, 'lb')+ logCcache(k, T, 'lb') - 2*n);
                    end
                end
            end
            
            function out = logCcache(n, k, UBorLB)
                if n < k || k < 0
                    out = 0;
                else
                    switch lower(UBorLB)
                        case 'ub'
                            out = logCub(n+1, k+1);
                        case 'lb'
                            out = logClb(n+1, k+1);
                    end
                end
            end
        end
        
        function out = bound(d, b)
            if ~isnan(boundCache(d))
                out = boundCache(d);
                return;
            end
            k = n;
            out = log2(b/sum(2.^( logCub(k+1, 1:floor(d) + 1) - n)));
            boundCache(d) = out;
            
        end
    end

%--------------------------------------------------------------------------
    function out = Approx()
        %Gaussian approximation (R = 1)
        %Vc = delta*(1-delta)*(log2((1-delta)/delta))^2;
        %Vs = p*(1-p)*(log2((1-p)/p))^2;
        lambda = log2((1-dbar)/dbar);
        
        out = dbar + sqrt((Vc + Vs)/n)/lambda*Qinv(e);
    end
%--------------------------------------------------------------------------
    function out = ARandom()
        %Random coding fair source
        out = fzero(@(d)Pe(d) - e, [delta 1]);
        
        
        function out = Pe(d)
            out = 0;
            for k = 0:ceil(n*d)
                out = out + binopdf(k, n, delta)*exp(- Csum(n, ceil(n*d - k), 'lb', accuracy)*2^(n*(R - 1)) );
            end
            
            out = out + 1 - binocdf(n*d, n, delta);
        end
    end

%--------------------------------------------------------------------------
%Utilities for faster computation
    function updatex0(x)
        if ~isnan(x) && n > 100
            x0 = x;
            nprev = n;
        end
    end


    function updatey0(y)
        if ~isnan(y)
            y0 = y;
        end
    end

end


