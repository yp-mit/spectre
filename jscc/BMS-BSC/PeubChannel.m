function out = PeubChannel(n, delta, L, accuracy)
%upper bound for the error probability of the binary symmetric channel 
%RCU bound. Holds for maximal error probability if M = 2^L where L is an
%integer

%n - block length (scalar)
%delta - channel crossover probability
%L = log M, where M is the number of messages
%accuracy - accuracy of computation (-1 for worse and 0 for better)


%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

if accuracy < 0 
    %normal approximation
    C = 1 - h(delta);
    Vc = delta*(1-delta)*(log2((1-delta)/delta))^2;
    out = Q( ( n*C - L + log2(n)/n )/sqrt(n*Vc) );  
    return;
end

out = 0;
for t = 0:n
    if bound(t) > 1
        out = out + 1 - binocdf(t - 1, n, delta); %Prob(X > t - 1) = Prob(X >= t)
        break;
    end
    out = out + binopdfbound(t, n, delta, 'ub')*bound(t);
end

    function out = bound(t)
        out = (2^(L-n) - 2^(-n))*Csum(n, t, 'ub', accuracy);
    end
end