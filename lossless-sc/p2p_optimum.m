function R = p2p_optimum(a, n_vals, eps)
% Compute the achievable sum rate (at the symmetrical rate point) based on the point-to-point optimal code 
% for the joint source distribution: [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)].

% Input:
% a: parameter to define the joint source distribution, should be in range [0.25, 1)
% n_vals: a vector containing the blocklengths
% eps: target error probability
% Output:
% R: the optimal rates for the blocklengths specified in n_vals.

    % Compute log binomial coefficients
    log_b = binomial_coeff(max(n_vals));

    R = zeros(length(n_vals),1);

    for i = 1:length(n_vals)
        n = n_vals(i);
        m = 0:n;  % Number of occurrence of the largest joint probability mass
        % log_Pr stores all the log probabilities in descending order
        log_Pr = m*log(a)+(n-m)*log(1/3*(1-a));  

        % Compute the optimal rate
        log_sum_Pr = log_b{n+1}(n+1)+log_Pr(n+1);
        k = n;
        while exp(log_sum_Pr) < 1 - eps   
            % Add up the probabilities of the strings until the total probability >= 1 - eps    
            k = k - 1;
            log_sum_Pr = my_logsumexp([log_sum_Pr log_b{n+1}(k+1)+(n-k)*log(3)+log_Pr(k+1)]);  
        end
        d = 1 - eps - (exp(log_sum_Pr) - exp(log_b{n+1}(k+1)+(n-k)*log(3)+log_Pr(k+1))); 
        M1 = my_logsumexp(log_b{n+1}(k+2:n+1)+(n-(k+1:n))*log(3));
        M2 = log(d) - log_Pr(k+1);
        R(i) = 1/n*my_logsumexp([M1 M2]);
    end
end

