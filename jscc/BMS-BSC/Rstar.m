function out = Rstar(p, delta, d, e, n, fun)
%finds the rate compatible with a given excess distortion at a given blocklength for coin
%flip source transmitted over a BSC
%p - source bias
%delta  - channel crossover probability
%d - excess distortion
%e - excess probability
%n - block length (scalar)
%fun - which function to use for calculation

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

%starting points - 'persistent' to make optimization faster
persistent x0;
persistent x0prev;
persistent nprev;
lenprev = 4;

if isempty(nprev)
    nprev = 0;
end

tol = 1e-12;
accuracy = 1; %accuracy: 1 - more accurate, 0 - faster

Capacity = 1 - h(delta);
Rd = h(p) - h(d);
Rlim = Capacity/Rd;
Vc = delta*(1-delta)*(log2((1-delta)/delta))^2;
Vs = p*(1-p)*(log2((1-p)/p))^2;
deltan = n - nprev; %distance between consecutive n's    
x0correction = 1/Rd*sqrt((Vc + Rlim*Vs))*Qinv(e)/2/n^(3/2)*deltan;
x0slack = 16*abs(x0correction);
y0slack = .2*Capacity;

%if length(x0prev) >= lenprev
%    x0slack = min(x0slack, 3*std(x0prev));
%end

if length(y0prev) >= lenprev
    y0slack = min(y0slack, 3*std(y0prev));
end

switch lower(fun)
    case 'approx'
        %gaussian approximation
        out = Approx();
    case 'clistdecoding'
        %Converse via list decoding
        if isempty(x0) || isempty(y0)
            %x0 = .8669; % p= 2/5 n = 541 %Rlim; %upper bound for rate
            x0 = .91; %p = 2/5 n = 940
            y0 = [0 1]; %for listdecodingfair
        end
        out = CListDecoding();
    case 'cgeneral'
        options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set', 'Display', 'off');
        if isempty(x0)
            %x0 = 1.6; %n = 350
        else
            x0 = x0 + x0correction;
        end
        
        out = Cgeneral();
    case 'approxseparate'
        %gaussian approximation
        options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set', 'Display', 'notify');
        out = ApproxSeparate();
    case 'aseparate'
        %Achievability via separate coding
        if isempty(x0) 
            %x0 = 1.3; % p =.11 n = 760
            %x0 = .84; %p = .4 n = 880
            %y0 = (Capacity + Rd*x0)/2;
        else       
            x0 = x0 + x0correction;
            y0 = (Capacity + Rd*x0)/2;
        end
        out = ASeparate();
    case 'ajoint'
        %Achievability via JSCC coding
        if isempty(x0) || isempty(y0)
        %    x0 = Rlim; %upper bound for rate
            y0 = 10; % gg
        end
        options = optimset('TolX',tol, 'MaxFunEvals', 500, 'Algorithm', 'active-set', 'Display', 'off');
        out = AJoint();   
    case 'arcu'
        %Achievability via RCU bound
        %if isempty(x0)
        %    x0 = Rlim; %upper bound for rate
        %end
        out = ARCU();
    otherwise
        disp('Error: unknown type.')
end


%--------------------------------------------------------------------------
    function out = Approx()
        %Gaussian approximation
        %Vc = delta*(1-delta)*(log2((1-delta)/delta))^2;
        if p == 1/2
            out = Rlim - sqrt(Vc/n)/Rd*Qinv(e);
        else
            %Vs = p*(1-p)*(log2((1-p)/p))^2;
            out = Rlim - 1/Rd*sqrt((Vc + Rlim*Vs)/n)*Qinv(e);
        end
    end

%--------------------------------------------------------------------------
    function out = CListDecodingFair()
        %converse via list decoding (fair source)
        
        %find rstar = max r: prob(r) <=1 - e
        rstar = fzero(@(r) binocdf(floor(r*n), n, delta) - 1 + e, y0);
        y0 = rstar;
        
        rstar = maxintless(@(r)binocdf(r, n, delta) - 1 + e, rstar*n);
        lambda = (1 - e - binocdf(rstar, n, delta))/binopdf(rstar +1, n, delta);
        
        kstar = fzero(@(r) bound(ceil(r*n), 0), [1/n x0]);
        %x0 = kstar + slack;
        kstar = maxintless(@(k)bound(k, accuracy), kstar*n);
        out = kstar / n;
        fprintf('CListDecoding: n = %i, R = %f\n', n, out);
        
        function out = bound(k, accuracy)
            out = log2( Csum(n, rstar, 'lb', accuracy) + lambda*C(n, rstar+1, 'lb') )...
                -  log2(Csum(k, floor(k*d), 'ub', accuracy ))  - n + k;
        end
    end

%--------------------------------------------------------------------------
    function out = CListDecoding()
        %Converse via meta-converse
        if p == 1/2
            out = CListDecodingFair();
            return;
        end
        
        if isempty(x0)
            %first pass
            startk = 1;
            endk = floor(Rlim*n);
        else
            startk = floor( max( 1, (x0-x0slack)*n ) );
            endk = floor( min( Rlim*n, (x0+x0slack)*n ) );
        end
        
        %precompute Binomial coeff
        sizeC = max(endk, n)+1;
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
        r = NaN(endk, n);
        lp = log2((1-p)/p);
        ldelta = log2((1 - delta)/delta);
        for k = 0:endk
            for j = 0:n
                r(k+1, j+1) = k*lp + j*ldelta;
            end
        end
        
        r = unique(sort(r(:)));
        
        boundCache = NaN(1, floor(Rlim*n));
        out = fzerointeger(@(k)-bound(k), startk, endk, 'min');
        out = out/n;
        updatex0(out);
        %fprintf('CListDecoding: n = %i, R = %f \n', n, out);
        
        
        function out = bound(k)
            if ~isnan(boundCache(k))
                out = boundCache(k);
                return;
            end
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
                        out = out + 2^(logCcache(n, t, 'lb')+logCcache(k, T, 'lb') - n);
                    elseif T*lp + t*ldelta - rstar < tol
                        out = out + lambda*2^(logCcache(n, t, 'lb')+ logCcache(k, T, 'lb') - n);
                    end
                end
            end
            out = log2(out/sum(2.^(logCub(k+1, 1:floor(k*d) + 1))));
            boundCache(k) = out;
            
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
    end

%--------------------------------------------------------------------------

    function out = Cgeneral()
        if isempty(x0)
            %first pass
            startk = 1;
            endk = floor(Rlim*n);
        else
            startk = floor( max( 1, (x0-x0slack)*n ) );
            endk = floor( min( Rlim*n, (x0+x0slack)*n ) );
        end
        
        taulb = tol;
        tauub = 10;
        tau0 = (taulb + tauub)/2; %tau (not updated because we need an array for each k)
        
        out = fzerointeger(@(x)-Pe(x/n) + e, startk, endk, 'min');
        out = out/n;
        updatex0(out);
        
        function out = Pe(R)
            shift = 5;
            [tau, f] = fmincon(@logPe, tau0,[],[],[],[],taulb, tauub, [], options);
            out = 2^(-f)-shift;
            %updatey0(tau);
            
            function out = logPe(tau)
                gamma = log2(n)/2*tau;
                k = round(R*n);
                out = 0;
                if p == 1/2
                    for c = 0:n
                        iXY = n*log2(2 - 2*delta) + c*log2(delta/(1 - delta));
                        jS = k*Rd;
                        if jS - iXY >= gamma
                            out = out + binopdfbound(c, n, delta, 'lb');
                        end
                    end
                else
                    for c = 0:n
                        for s = 0:k
                            iXY = n*log2(2 - 2*delta) + c*log2(delta/(1 - delta));
                            jS = -s*log2(p) - (k - s)*log2(1 - p) - k*h(d);
                            if jS - iXY >= gamma
                                out = out + binopdfbound(c, n, delta, 'lb')*binopdfbound(s, k, p, 'lb');
                            end
                        end
                    end
                end
                out = out - 2^(-gamma);
                out = - log2(shift + out);
            end
        end
    end

%--------------------------------------------------------------------------
    function out = AJoint()
        if isempty(x0)
            %first pass
            startk = 1;
            endk = floor(Rlim*n);
        else
            startk = floor( max( 1, (x0-x0slack)*n ) );
            endk = floor( min( Rlim*n, (x0+x0slack)*n ) );
        end
        
        %precompute Binomial coeff
        sizeC = max(endk, n)+1;
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
        
        
        boundCache = NaN(1, floor(Rlim*n));
        out = fzerointeger(@(k)bound(k) - e, startk, endk, 'max');
        out = out/n;
        updatex0(out);
%       fprintf('AJoint: n = %i, R = %f \n', n, out);
        
        function out = bound(k)
            if ~isnan(boundCache(k))
                out = boundCache(k);
                return;
            end
            
            %optimize with respect to r and g; y0 = [r gg]
            [gg, out] = fmincon(@goal,y0,[],[],[],[], 0, Inf, [], options);
            
            boundCache(k) = out;
            
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
                    T = ceil((p-d)/(1-2*d)*k);
                    j0 = ceil( (T + s - k*d)/2 );
                    out = 0;
                    for j = j0:s
                        out = out + 2^(logCcache(T, j, 'lb') + logCcache(k-T, s-j, 'lb') - logCcache(k, s, 'ub'));
                    end
                end
                
                function out = Pdballfair(k, d)
                    out = 2^(-k)*Csum(k, floor(k*d), 'lb', accuracy);
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
        
        startk = max( ceil(n*(x0 - x0slack/2)), 1);
        endk = min( ceil(n*(x0 + x0slack)), floor(Rlim*n) );
           
        startL = max( ceil((y0 - y0slack)*n), 1);
        endL = min( ceil((y0 + y0slack)*n), n );
        
        if isempty(x0) || n <= 100
            %first pass
            startk = 1;
            endk = floor(Rlim*n);
        end
        
        if isempty(y0) || n <= 100
            startL = 1;
            endL = n;
        end
        
        %precompute channel error probability
        
        %At first find interesting region
        stepL = floor((endL - startL)/8);
        L = startL:stepL:endL;
        PeCh = Inf(1, numel(L));
        for i = 1:numel(L) %2^(L/step) - the number of messages passed from the source enc to the channel enc
            PeCh(i) = PeubChannel(n, delta, L(i), accuracy);
        end
        
        indL = find(PeCh <=e);
        
        if isempty(indL)
            out = NaN;
            return;
        end
        
        %refine
        startL = L(indL(1)) - stepL; 
        endL = L(indL(numel(indL)))+ stepL;
        stepL = 1;
        L = startL:stepL:endL;
        PeCh = Inf(1, numel(L));
        for i = 1:numel(L) %2^(L/step) - the number of messages passed from the source enc to the channel enc
            PeCh(i) = PeubChannel(n, delta, L(i), accuracy);
        end
        
        indL = find(PeCh>= e/100 & PeCh <=e);
        
        if isempty(indL)
            out = NaN;
            return;
        end
        
        kstar  = fzerointeger(@(k)Peubseparate(k, p, d, PeCh(indL), L(indL), accuracy) - e, startk, endk, 'max');
        out = kstar/n;        
        updatex0(out);
    end

%--------------------------------------------------------------------------

    function out = ApproxSeparate()
        %Gaussian approximation of separate source-channel coding
        %Vc = delta*(1-delta)*(log2((1-delta)/delta))^2;
        if p == 1/2
            out = Rlim - sqrt(Vc/n)/Rd*Qinv(e);
        else
            %Vs = p*(1-p)*(log2((1-p)/p))^2;
            [~, opt] = fmincon(@goal,[e/2, e/2],[1 1], e,[],[], [0 0], [e e], [], options);
            out = Rlim - 1/sqrt(n)/Rd*opt;
        end
        
        function out = goal(arg)
            out = sqrt(Vs*Rlim)*Qinv(arg(1)) + sqrt(Vc)*Qinv(arg(2));
        end
    end


%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%% LOSSLESS achievability bounds %%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
    function out = ARCU()
        %Achievability (RCU) for d = 0
        
        if isempty(x0)
            %first pass
            startk = 1;
            endk = floor(Rlim*n);
        else
            startk = floor((x0-x0slack)*n );
            endk = floor((x0+x0slack)*n );
        end
        
        %precompute Binomial coeff
        sizeC = max(endk, n)+1;
        Cub = NaN(sizeC, sizeC);
        for k = 0:sizeC
            for j = 0:sizeC
                Cub(k+1, j+1) = C(k, j, 'ub');
            end
        end
        
        [out Peub] = fzerointeger(@(k)bound(k) - e, startk, endk, 'max');
        out = out/n;
        updatex0(out);
        fprintf('ARCU: n = %i, R = %f, Peub = %f\n', n, out, Peub);
        
        function out = bound(k)
            out = 0;
            for t = 0:n
                for T = 0:k
                    err = 0;
                    for j = 0:k
                        arg = min(n, floor(t + (T-j)*(log2(1-p) - log2(p))/(log2(1-delta) - log2(delta))));
                        err =  err + Ccache(k, j)/2^n * sum( Cub(n+1, 1:arg +1 ) );
                        if err >= 1
                            err = 1;
                            break;
                        end
                    end
                    out = out + Ccache(n, t)*delta^t*(1-delta)^(n-t)*Ccache(k, T)*p^T*(1-p)^(k-T)*err;
                end
            end
            function out = Ccache(n, k)
                if n < k || k < 0
                    out = 0;
                else
                    out = Cub(n+1, k+1);
                end
            end
            
        end
        
    end

%--------------------------------------------------------------------------
%Utilities for faster computation
    function updatex0(x)
        if ~isnan(x)
            x0 = x;
            nprev = n;
            
            if length(x0prev) < lenprev
                x0prev = [x0prev x0];
            else
                x0prev(1:lenprev - 1) = x0prev(2:lenprev);
                x0prev(lenprev) = x0;
            end
        end
    end

end
