function R = p2p_rcu(a, n_vals, eps)
% Compute the achievable sum rate (at the symmetrical rate point) from the point-to-point RCU bound 
% for the joint source distribution: [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)].

% Input:
% a: parameter to define the joint source distribution, should be in range [0.25, 1)
% n_vals: a vector containing the blocklengths
% eps: target error probability
% Output:
% R: the achievable rates for the blocklengths specified in n_vals.

    p = [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)];
    H = sum(sum(p.*log(1./p)));  % Source joint entropy

    % Compute log binomial coefficients
    log_b = binomial_coeff(max(n_vals));

    R = zeros(length(n_vals),1);

    for i = 1:length(n_vals)
        n = n_vals(i);
        m = 0:n;  % Number of occurrence of the largest joint probability mass
        % log_Pr stores all the log probabilities in descending order
        log_Pr = m*log(a)+(n-m)*log(1/3*(1-a));  

        % [x, y]: initial range for the bisection algorithm
        x = H;
        y = H + 2;
        %%%%%%%%%% Use this part to check if the range is valid for bisection algorithm %%%%%%%%%%
        c = x;
        log_indicator = zeros(1,n+1);   
        for j = 0:n   
            log_indicator(j+1) = min(0, my_logsumexp(log_b{n+1}(j+1:n+1)+(n-(j:n))*log(3))-n*c);
        end
        err = exp(my_logsumexp(log_Pr + log_b{n+1} + (n-m)*log(3) + log_indicator));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Start of bisection iteration
        iter = 0;
        while abs(x - y) > 0.00001 && iter <= 25
            c = (x+y)/2;
            %%%%%%%%%%%%%%%% body of iteration %%%%%%%%%%%%%%%%
            log_indicator = zeros(1,n+1);   
            for j = 0:n   
                log_indicator(j+1) = min(0, my_logsumexp(log_b{n+1}(j+1:n+1)+(n-(j:n))*log(3))-n*c);
            end
            err = exp(my_logsumexp(log_Pr + log_b{n+1} + (n-m)*log(3) + log_indicator));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if err - eps < 0
                y = c;
            else
                x = c;
            end
            iter = iter + 1;
        end
        R(i) = c;
    end
end

