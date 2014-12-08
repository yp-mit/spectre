function f = kappa_inf(tau, P)
%  This function computes approximation to kappa in the limit of n->inf
%
x0 = norminv( (tau+1)/2, 0, 1);
VP = 2 * (1 + 2*P);
VQ = 2 * (1 + P)^2;
f = 2*normcdf(sqrt(VP/VQ) * x0, 0, 1) - 1;

