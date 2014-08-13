function f = kappa_inf(tau, A)
%  This function computes approximation to kappa in the limit of n->inf
%
x0 = norminv( (tau+1)/2, 0, 1);
VP = 2 * (1 + 2*A^2);
VQ = 2 * (1 + A^2)^2;
f = 2*normcdf(sqrt(VP/VQ) * x0, 0, 1) - 1;

