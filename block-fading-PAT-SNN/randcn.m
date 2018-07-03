function z = randcn(m,n)
%CRANDN returns samples from a complex normal random variable with unit variance
%   Z = RANDCN(M,N) returns an m-by-n complex normal random matrix with column
%   vectors with identity covariance.
%   Z = RANDCN returns one realisation only.

if nargin == 0
    m = 1;
    n = 1;
end

if nargin == 1
    n = m;
end

z = sqrt(0.5) * ( randn(m,n) + 1i*randn(m,n) ); 

end
    