function i_s = idgallagerqpsknn(x,y,hest,s,rho)
%IDGALLAGERNN returns the information density with Gallager exponent s
%   IS = IDGALLAGERQPSKNN(X,Y,HEST,S) returns the information density 
%       i_s(x,y,hest) = log( q(x,y,hest)^s / E_X[q(X,y,hest)^s] )
%   where q(x,y,hest) is the channel believed true by the nearest neighbor (NN) 
%   decoder.
%
%   N.B. All parameters are assumed scalars.
%       

% Compute q(x,y,hest) up to a normalisation constant
q = exp( - abs( y-hest*x )^2 );

% Compute E_X[q(X,y,hest)^s] up to a normalisation constant
qs_avg = 0;
for k=1:4
    x_k = sqrt(rho) * exp(1i*pi/2 * (k-1));
    qs_avg = qs_avg + exp( - s * abs( y-hest*x_k )^2 );
end
qs_avg = qs_avg/4;
%qs_avg = mean( exp( - s * abs( y-hest*sqrt(rho)*exp(1i*pi/2*(0:3)) ).^2 ) );

i_s = s * log( q ) - log( qs_avg ); 