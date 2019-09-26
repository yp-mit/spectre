function R = sw_han_converse(a, n_vals, eps)
% Compute the converse sum rate (at the symmetrical rate point) from Han's SW converse 
% for the joint source distribution: [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)].

% Input:
% a: parameter to define the joint source distribution, should be in range [0.25, 1)
% n_vals: a vector containing the blocklengths
% eps: target error probability
% Output:
% R: the rates for the blocklengths specified in n_vals.

    p = [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)];                                 % Stores joint probabilities
    H = sum(sum(p.*log(1./p)));                                             % Joint entropy
    b = [a/(a+1/3*(1-a)) 1/3*(1-a)/(2*1/3*(1-a)) 1/3*(1-a)/(a+1/3*(1-a))];  % Stores conditional probabilities

    % Compute log binomial coefficients
    log_b = binomial_coeff(max(n_vals));

    R = zeros(length(n_vals),1);

    for i = 1:length(n_vals)
        n = n_vals(i);
        m = 0:n;
        log_Pr = m*log(a)+(n-m)*log(1/3*(1-a));
        
        % Search the range of gamma specified in r_vals
        % The range of r_vals could be set smaller for large n to save computation time.
        if n <= 500
            r_vals = 0.0001:0.00005:2;
        else 
            r_vals = 0.0001:0.00005:0.1;
        end
        
        R_r = zeros(1,length(r_vals));
        for ri = 1:length(r_vals)
            r = r_vals(ri);
            x = 0;
            y = H + 0.1;
            
            %%%%%%%%%% Use this part to check if the range is valid for bisection algorithm %%%%%%%%%%
            c = x;
            m0 = max(0, floor(n*(log(1/3*(1-a))+c+r)/log(1/3*(1-a)/a)));
            if m0 >= n
                err = 1 - 3*exp(-n*r);
            else
                log_err_m1 = my_logsumexp(log_Pr(1:m0+1)+log_b{n+1}(1:m0+1)+(n-(0:m0))*log(3));
                log_err_m2 = [];
                for m = m0+1:n
                    B_m = ceil((m*log(b(1))+(n-m)*log(b(2))+n*(c/2+r))/(log(b(2)/b(3))));  
                    if B_m <= 0
                        C_m = log_Pr(m+1) + log_b{n+1}(m+1) + (n-m)*log(3);
                        log_err_m2 = [log_err_m2 C_m];
                    elseif B_m <= n-m
                        C_m1 = my_logsumexp(log_b{n-m+1}(B_m+1:n-m+1)+(n-m-(B_m:n-m))*log(2));               
                        log_err_k = [];
                        for k = 0:B_m-1
                            if B_m <= n-m-k
                                D_m = log_b{n-m+1}(k+1) + my_logsumexp(log_b{n-m-k+1}(B_m+1:n-m-k+1));
                                log_err_k = [log_err_k D_m];
                            end
                        end
                        C_m2 = my_logsumexp(log_err_k);
                        C_m = log_Pr(m+1) + log_b{n+1}(m+1) + my_logsumexp([C_m1 C_m2]);        
                        log_err_m2 = [log_err_m2 C_m];
                    end          
                end
                log_err = my_logsumexp([log_err_m1 log_err_m2]);      
                err = exp(log_err) - 3*exp(-n*r);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Start of bisection iteration
            iter = 0;
            while abs(x - y) >= 0.00001 && iter <= 30
                c = (x+y)/2;
                %%%%%%%%%%%%%%%%%%% body of iteration %%%%%%%%%%%%%%%%%%%%
                m0 = max(0, floor(n*(log(1/3*(1-a))+c+r)/log(1/3*(1-a)/a)));
                if m0 >= n
                    err = 1 - 3*exp(-n*r);
                else
                    log_err_m1 = my_logsumexp(log_Pr(1:m0+1)+log_b{n+1}(1:m0+1)+(n-(0:m0))*log(3));
                    log_err_m2 = [];
                    for m = m0+1:n
                        B_m = ceil((m*log(b(1))+(n-m)*log(b(2))+n*(c/2+r))/(log(b(2)/b(3))));  
                        if B_m <= 0
                            C_m = log_Pr(m+1) + log_b{n+1}(m+1) + (n-m)*log(3);
                            log_err_m2 = [log_err_m2 C_m];
                        elseif B_m <= n-m
                            C_m1 = my_logsumexp(log_b{n-m+1}(B_m+1:n-m+1)+(n-m-(B_m:n-m))*log(2));               
                            log_err_k = [];
                            for k = 0:B_m-1
                                if B_m <= n-m-k
                                    D_m = log_b{n-m+1}(k+1) + my_logsumexp(log_b{n-m-k+1}(B_m+1:n-m-k+1));
                                    log_err_k = [log_err_k D_m];
                                end
                            end
                            C_m2 = my_logsumexp(log_err_k);
                            C_m = log_Pr(m+1) + log_b{n+1}(m+1) + my_logsumexp([C_m1 C_m2]);        
                            log_err_m2 = [log_err_m2 C_m];
                        end          
                    end
                    log_err = my_logsumexp([log_err_m1 log_err_m2]);      
                    err = exp(log_err) - 3*exp(-n*r);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if err - eps < 0
                    y = c;
                else
                    x = c;
                end
                iter = iter + 1;
            end
            R_r(ri) = c;
        end
        % Optimize with respect to gamma
        R(i) = max(R_r);
    end
end




