function x = qpsksample(m,n)
%QPSKSAMPLE returns QPSK symbols on the unit circle starting with (1,0)
%   X = QPSKSAMPLE(M,N) returns an m-by-n matrix of iid symbols 
%   sampled uniformly at random from a QPSK constellation.
%   X = QPSKSAMPLE returns one symbol only

if nargin == 0
    m = 1;
    n = 1;
end

if nargin == 1
    n = 1;
end

x = exp(1i*pi/2 * randi([1,4],m,n) ); 

end