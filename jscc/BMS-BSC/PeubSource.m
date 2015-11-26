function out = PeubSource(k, p, d, L, accuracy)
%upper bound for the error probability of the binary source
%achievability via sphere covering
%k - source blocklength
%p - source bias
%d - distortion threshold
%L = log M, where M is the number of representation points
%accuracy = accuracy of Csum: 1 - more accurate, 0 - faster

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

Kbig = 80;

if accuracy < 0
    %normal approximation
    Rd = h(p) - h(d);
    Vs = p*(1-p)*(log2((1-p)/p))^2;
    if p ~= 1/2 && k >= Kbig
        out = Q( (L - k*Rd - log2(k))/sqrt(k*Vs) );
    else
        out = PeubSource(L, k, p, d);
    end
    return;
end

if d == 0
    %lossless
    out = PeubSourceLossless(L, k, p);
elseif p == 1/2
    %FAIR coin flip source
    out = PeubSourceFair(L, k, d);
else
    %arbitrary p source
    out = PeubSourceBiased(L, k, p, d);
end
%        fprintf('SphereCoveringA: L = %i, k = %i, Peub = %f\n', L, k, out);

    function out = PeubSourceLossless(L, k, p)
        csum = 0;
        kstar = k;
        for i = 0:k
            csum = csum + C(k, i, 'ub');
            if log2(csum) > L
                kstar = i - 1;
                break;
            end
        end
        out =  1 - binocdf(kstar,k,p);
    end

    function out = PeubSourceFair(L, k, d)
        if k <= Kbig
            out =  2^L *log(1 - Csum(k, floor(k*d), 'lb', accuracy) /(2^k));
        else
            out = -2^(L-k)*Csum(k, floor(k*d), 'lb', accuracy);
        end
        out = exp(out);
    end

    function out = PeubSourceBiased(L, k, p, d)
        q = max ( 0, (p-d)/(1-2*d));
        if q == 0
            q = p;
        end
        out = 0;
        for t = 0:k
            out = out + binopdfbound(t,k,p, 'ub')*e(k,t);
            if isnan(out)
                break;
            end
        end
        
        function out = e(k,t)
            out = 0;
            if 1
                for T = 0:k
                    out = out + binopdfbound(T, k, q, 'ub')*Pwithin(k, t, T);
                end
            else
                T = round(q*k); %#ones
                out = Pwithin(k, t, T);
            end
            
            if k <= Kbig
                out = 2^L*log( 1 - out);
                out = exp(out);
            else
                out = exp(-2^L*out);
            end
            
        end
        
        function out = Pwithin(k,t, T)
            %y of type T is within a given x of type t - almost exact lower bound
            t0 = max ( 0, ceil( (t+T - k*d)/2) ); %# ones  in 11111 zone (T ones in codewords)
            
            out = logC(T, t0, 'lb') + logC(k-T, t - t0, 'lb') - logC(k, t, 'ub');
            out = 2^out;
        end
        
    end
end