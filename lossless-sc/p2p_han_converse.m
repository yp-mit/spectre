function R = p2p_han_converse(a, n_vals, eps)
% Compute the converse sum rate (at the symmetrical rate point) from Han's point-to-point converse 
% for the joint source distribution: [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)].

% Input:
% a: parameter to define the joint source distribution, should be in range [0.25, 1)
% n_vals: a vector containing the blocklengths
% eps: target error probability
% Output:
% R: the rates for the blocklengths specified in n_vals

    p = [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)];      % Stores joint probabilities
    H = sum(sum(p.*log(1./p)));                                             

    % Compute log binomial coefficients
    log_b = binomial_coeff(max(n_vals));

    R = zeros(length(n_vals),1);

    for i = 1:length(n_vals)
        n = n_vals(i);
        m = 0:n;
        log_Pr = m*log(a)+(n-m)*log(1/3*(1-a));  
        
        % Search the range of gamma specified in gamma_vals
        % The range of gamma_vals could be set smaller for large n to save computation time.
        gamma_vals = 0.0001:0.00005:0.5;
        
        R_gamma = zeros(1,length(gamma_vals));
        for ri = 1:length(gamma_vals)
            r = gamma_vals(ri);
            % [x, y]: initial range for the bisection algorithm
            x = 0;
            y = H + 2; 
            %%%%%%%%%% Use this part to check if the range is valid for bisection algorithm %%%%%%%%%%
            c = x;      
            m0 = max(0, floor(n*(log(1/3*(1-a))+c+r)/log(1/3*(1-a)/a)));
            if m0 >= n
                err = 1 - exp(-n*r);
            else
                log_err = my_logsumexp(log_Pr(1:m0+1)+log_b{n+1}(1:m0+1)+(n-(0:m0))*log(3));    
                err = exp(log_err) - exp(-n*r);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Start of bisection iteration
            iter = 0;
            while abs(x - y) >= 0.00001 && iter <= 30
                c = (x+y)/2;
                %%%%%%%%%%%%%%%%%%% body of iteration %%%%%%%%%%%%%%%%%%%%
                m0 = max(0, floor(n*(log(1/3*(1-a))+c+r)/log(1/3*(1-a)/a)));
                if m0 >= n
                    err = 1 - exp(-n*r);
                else
                    log_err = my_logsumexp(log_Pr(1:m0+1)+log_b{n+1}(1:m0+1)+(n-(0:m0))*log(3));    
                    err = exp(log_err) - exp(-n*r);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if err - eps < 0
                    y = c;
                else
                    x = c;
                end
                iter = iter + 1;
            end
           R_gamma(ri) = c; 
        end
        % Optimize with respect to gamma
        R(i) = max(R_gamma);
    end
end


