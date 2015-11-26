function out = PeubChannel(P, n, L)
%upper bound for the error probability of the AWGN channel
%P - channel SNR
%n - block length (scalar)
%L - length of messages transmitter between the source encoder and the channel encoder   

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%


%dP/dQ is bounded by a constant. Find that constant.
[~, gamma2] = fminbnd(@(x)-ncx2pdf(n*x, n, n*P)./chi2pdf(n*x/(1+P), n)*(1+P), 0, 10);
gamma2 = -gamma2;

A = n/2*log(1+P) + n/2 - L - log(gamma2);


Finner = @(v) ncx2pdf(n*v, n, n/P).*n.*exp(P/2/(1+P)*n*(v - thresv()));
tail = 1 - ncx2cdf(thresv()*n,n,n/P);
out = quad(Finner, 0, thresv())...                      %exponent < 1
    + tail;                                             %exponent = 1

    function out = thresv()
        %threshold for v as a function of v0
        out = max( 0, 2*(1+P)/P/n*A );
    end
end