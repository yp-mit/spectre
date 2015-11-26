function [out, Lstar] = Peubseparate(k, n, p, delta, d, startL, endL, accuracy, options)
%Excess-distortion probability for the transmission of a BMS
%over a BSC as achieved by a separated scheme - optimized with respect to the number of
%messages allowable between the channel encoder and the channel decoder (2^L)
%k - source blocklength
%n - channel blocklength (scalar)
%p - source bias
%delta - channel crossover probability
%d - distortion threshold
%startL - lowest L considered
%endL - highest L considered
%accuracy - -1 for worse and 0 for better
%options - optimization options

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

[Lstar, out] = fminbnd(@(L)PeubSource(k, p, d, L, accuracy) + PeubChannel(n, delta, L, accuracy),startL, endL, options);
tmp = round(2^Lstar);
Lstar = log2(tmp);

end






