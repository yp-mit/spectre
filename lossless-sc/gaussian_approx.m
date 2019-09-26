function R = gaussian_approx(a, n_vals, eps)
% Compute the third-order Gaussian approximation of the sum rate (at the symmetrical rate point)
% for the joint source distribution: [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)].

% Input:
% a: parameter to define the joint source distribution, should be in range [0.25, 1)
% n_vals: a vector containing the blocklengths
% eps: target error probability
% Output:
% R: the third-order Gaussian approximation for the blocklengths specified in n_vals.

    p = [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)];
    H = sum(sum(p.*log(1./p)));
    V = sum(sum(p.*(log(1./p)-H).^2));
    R = H + sqrt(V./n_vals)*qfuncinv(eps) - log(n_vals)./(2*n_vals);
end

