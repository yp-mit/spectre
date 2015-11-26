function [out, bestL] = Peubseparate(k, d, PeCh, e)
%Excess-distortion probability for the transmission of a Gaussian source
%over a Gaussian channel as achieved by a separated scheme - optimized with respect to the number of
%messages allowable between the channel encoder and the channel decoder (2^L)

%PeCh - array of channel error probabilities
%k - source blocklength
%d - target distortion
%e - excess distortion probability

%outputs: 
%out - the excess distortion probability
%bestL - the optimal length of messages exchanged between the encoder and
%the decoder

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

Kbig = 20;
tol = 1e-3;

out = Inf;
bestL = NaN;
for L = find(PeCh <= e)
    cur = PeubSource(L, k, d) + PeCh(L);
    if out > cur
        out = cur;
        bestL = L;
    end
end
%fprintf('ASeparate: k = %i, Peub = %f\n', k, out);

    function out = PeubSource(L, k, d)
        %achievability via sphere covering
        
        r0 = sqrt( max( 0, 1 - d) );
        a = r0 - sqrt(d);
        b = r0 + sqrt(d);
        
        out = chi2cdf((a^2 + tol)*k, k) + 1 - chi2cdf((b^2 - tol)*k, k)...  %nontypical source realization
            + quad(@Pes, a^2 + tol, b^2- tol);                              %source error
        
        function out = logPdball(x)
            %x - distance squared from the origin (x*k central chi square k)
            cosa = (r0^2 + x - d)./(2*r0.*sqrt(x));
            sina = (1 - cosa.^2).^(1/2);
            out = loggamma(k/2+1, 'LB') - loggamma((k-1)/2 +1, 'UB')- .5*log(pi)-log(k) + (k-1)*log(sina);
        end
        
        function out = Pes(x)
            %source error conditioned on the distance squared from the origin (x*k central chi square k)
            if k <  Kbig
                out = (1 - exp(logPdball(x))).^exp(L);
            else
                out = exp( - exp( logPdball(x) + L));
            end
            out = out.*chi2pdf(k*x, k).*k;
        end
    end
end



