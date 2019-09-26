function R = sw_rcu(a, n_vals, eps)
% Compute the achievable sum rate (at the symmetrical rate point) from the SW RCU bound 
% for the joint source distribution: [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)].

% Input:
% a: parameter to define the joint source distribution, should be in range [0.25, 1)
% n_vals: a vector containing the blocklengths
% eps: target error probability
% Output:
% R: the achievable rates for the blocklengths specified in n_vals.

    p = [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)];
    H = sum(sum(p.*log(1./p)));

    % Compute log binomial coefficients
    log_b = binomial_coeff(max(n_vals));

    R = zeros(length(n_vals),1);

    for i = 1:length(n_vals)
        n = n_vals(i);
        m = 0:n;
        % log_Pr stores all the log joint probabilities in ascending order
        log_Pr = m*log(a)+(n-m)*log(1/3*(1-a));  

        % Calculate and store the coefficients that will be used iteratively later 
        A_m = {};
        A_kj = {};
        for m = 0:n
            A_tempt = log_b{n+1}+(n:-1:0)*log(3);
            A_m{m+1} = my_logsumexp(A_tempt(m+1:n+1));
            for k = 0:n-m
                A_kj{m+1}(k+1) = my_logsumexp(log_b{m+k+1}(m+1:m+k+1))+(n-m-k)*log(2);
            end
        end

        % [x, y]: initial range for the bisection algorithm
        % For n small or eps small, need to adjust the right boundary of this initial interval. 
        x = H;
        y = H + 0.5;
        %%%%%%%%%% Use this part to check if the range is valid for bisection algorithm %%%%%%%%%%
        c = x;
        log_err_m = zeros(1,n+1);
        for m = 0:n
            log_err_k = zeros(1,n-m+1);
            for k = 0:n-m
                log_err_j = zeros(1,n-m-k+1); 
                for j = 0:n-m-k
                    log_err_j(j+1) = log_b{n-m-k+1}(j+1) + ...
                        min(0, my_logsumexp([A_kj{m+1}(j+1)-n*c/2 A_kj{m+1}(k+1)-n*c/2 A_m{m+1}-n*c]));
                end
                log_err_k(k+1) = log_b{n-m+1}(k+1) + my_logsumexp(log_err_j);
            end
            log_err_m(m+1) = log_Pr(m+1) + log_b{n+1}(m+1) + my_logsumexp(log_err_k);
        end
        err = exp(my_logsumexp(log_err_m)); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Start of bisection iteration
        iter = 0;
        while abs(x - y) > 0.00001 && iter <= 30
            c = (x+y)/2;     
            %%%%%%%%%%%%%%%% body of iteration %%%%%%%%%%%%%%%%
            log_err_m = zeros(1,n+1);
            for m = 0:n
                log_err_k = zeros(1,n-m+1);
                for k = 0:n-m
                    log_err_j = zeros(1,n-m-k+1); 
                    for j = 0:n-m-k
                        log_err_j(j+1) = log_b{n-m-k+1}(j+1) + ...
                            min(0, my_logsumexp([A_kj{m+1}(j+1)-n*c/2 A_kj{m+1}(k+1)-n*c/2 A_m{m+1}-n*c]));
                    end
                    log_err_k(k+1) = log_b{n-m+1}(k+1) + my_logsumexp(log_err_j);
                end
                log_err_m(m+1) = log_Pr(m+1) + log_b{n+1}(m+1) + my_logsumexp(log_err_k);
            end
            err = exp(my_logsumexp(log_err_m)); 
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




