function [R, log_c_binary] = p2p_ht_converse(a, n_vals, eps)
% Compute the converse sum rate (at the symmetrical rate point) from the binary hypothesis testing
% converse for the joint source distribution: [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)]. 
% This is a helper function called by sw_ht_converse to compute initial values for 
% the parameters in computing the Slepian-Wolf composite hypothesis testing converse. 

% Input:
% a: parameter to define the joint source distribution, should be in range [0.25, 1)
% n_vals: a vector containing the blocklengths
% eps: target error probability
% Output:
% R: converse sum rate, which is equivalent to the point-to-point optimum sum rate computed by p2p_optimum.
% log_c_binary: the initial values of the parameters for the blocklengths specified in n_vals

    % Compute log binomial coefficients
    log_b = binomial_coeff(max(n_vals));

    R = zeros(length(n_vals),1);
    log_c_binary = zeros(length(n_vals),1);

    for i = 1:length(n_vals)
        n = n_vals(i);
        m = 0:n; 
        log_Pr = m*log(a)+(n-m)*log(1/3*(1-a)); 

        % Compute the optimal rate
        k = n;
        log_sum_Pr = log_b{n+1}(k+1)+log_Pr(k+1);
        while exp(log_sum_Pr) < 1 - eps   
            k = k - 1;
            log_sum_Pr = my_logsumexp([log_sum_Pr log_b{n+1}(k+1)+(n-k)*log(3)+log_Pr(k+1)]);  % Add up the probabilities of the strings until the total probability >= 1 - eps    
        end
        log_c_binary(i) = log_Pr(k+1)/n;
        d = 1 - eps - (exp(log_sum_Pr) - exp(log_b{n+1}(k+1)+(n-k)*log(3)+log_Pr(k+1))); 
        log_lambda = log(d) - (log_Pr(k+1) + (n-k)*log(3));    
        U1 = my_logsumexp(log_b{n+1}(k+2:n+1)+(n-(k+1:n))*log(3));
        U2 = log_lambda + (n-k)*log(3);
        R(i) = 1/n*my_logsumexp([U1 U2]);
    end
end

