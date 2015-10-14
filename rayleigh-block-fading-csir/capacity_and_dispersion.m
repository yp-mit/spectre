function [C,V] = capacity_and_dispersion(n_t,n_r,T,P)

% This function numerically computes the capacity and dispersion of the 
% complex MIMO block fading channel with parameters
%
%   n_t - number of transmit antennas
%   n_r - number of reciever antennas
%   T   - Coherence time 
%   P   - Power constraint
%
% Note that one term in the dispersion is computed via monte carlo because
% it's an integral over a joint distribution of the eigenvalues of a
% Wishart W(I_{n_r},n_t) distribution.

lower = 0;
upper = Inf;

t = min(n_r,n_t);
r = max(n_r,n_t);

cap = integral(@(x) capacity_integrand(x,n_t,n_r,P,r,t),lower,upper);
C = min(n_t,n_r)*cap;

var_cap_term = monte_carlo_var_cap(n_t,n_r,P,T);
awgn_disp = integral(@(x) awgn_disp_integrand(x,n_t,n_r,P,T,r,t),lower,upper);
xi1 = integral(@(x) x1_integrand(x,n_t,n_r,P,T,r,t),lower,upper);
xi2 = integral(@(x) x2_integrand(x,n_t,n_r,P,T,r,t),lower,upper);

% For MISO, the expression for dispersion is different
if (n_r == 1) 
    if (n_t >= 2)
        disp('Warning: for MISO (n_r = 1), the calculated dispersion is an upper bound, can be improved with orthogonal design inputs')
        final_term = t/n_t*xi2^2;
    end
else
    final_term = t/n_t*xi2^2;
end

V = T*var_cap_term + (t*(1-awgn_disp) + t*(P/n_t)^2*(xi1 - final_term))*log2(exp(1))^2;
end

% E[1 - 1/(1 + P/n_t*lambda^2)^2]
function val = awgn_disp_integrand(lambda,n_t,n_r,P,T,r,t)
    % lambda_norm is the eigenvalues of 1/n_t*XX^T, where X is iid N(0,1)
    lambda_norm = 1/n_t*lambda;
    val =(1./(1+P/n_t*lambda_norm).^2).*lambda_marginal(lambda,r,t);
end

% E[(lambda^2 / (1 + P/n_t*lambda^2))^2]
function val = x1_integrand(lambda,n_t,n_r,P,T,r,t)
    lambda_norm = 1/n_t*lambda;
    val = ((lambda_norm./(1+P/n_t*lambda_norm)).^2).*lambda_marginal(lambda,r,t);
end

% E[lambda / (1 + P/n_t*\lambda^2)]^2
function val = x2_integrand(lambda,n_t,n_r,P,T,r,t)
    lambda_norm = 1/n_t*lambda;
    val = (lambda_norm./(1+P/n_t*lambda_norm)).*lambda_marginal(lambda,r,t);
end

% Var( \sum_{i=1}^{n_min} log(1 + P\n_t HH^T) )
function val = monte_carlo_var_cap(n_t,n_r,P,T)
    samples = 10^5;
    vals = zeros(1,samples);
    for i = 1:samples
        h = 1/sqrt(n_t)*sqrt(.5)*(randn(n_r,n_t) + 1i*randn(n_r,n_t));
        vals(i) = real(log2(det(eye(n_r) + P/n_t*h*h')));
    end
    val = (samples-1)/samples*var(vals);
end

% E[ log(1 + P/n_t HH^T) ]
function val = capacity_integrand(lambda,n_t,n_r,P,r,t)
    lambda_norm = 1/n_t*lambda;
    val = log2(1+P/n_t.*lambda_norm).*lambda_marginal(lambda,r,t);
end

% Computes the marginal distribution of an eigenvalue of a random matrix
% with Wishart(n_t,I_{n_r}) distribution
function p_lambda = lambda_marginal(lambda,r,t)
    to_sum = [];
    p_lambda = zeros(1,length(lambda));
    for i=1:length(lambda)
        for k=0:(t-1)
            to_sum = [to_sum factorial(k)/(factorial(k+r-t))*laguerreL(k,r-t,lambda(i))^2*lambda(i)^(r-t)*exp(-lambda(i))];
        end
        p_lambda(i) = 1/t*sum(to_sum);
        to_sum = [];
    end
end