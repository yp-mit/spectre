function R = sw_ht_converse(a, n_vals, eps)
% Compute the converse sum rate (at the symmetrical rate point) from Han's point-to-point converse 
% for the joint source distribution: [a 1/3*(1-a); 1/3*(1-a) 1/3*(1-a)].

% Input:
% a: parameter to define the joint source distribution, should be in range [0.25, 1)
% n_vals: a vector containing the blocklengths
% eps: target error probability
% Output:
% R: lower bounds of the rates for the blocklengths specified in n_vals
% Other:
% One can also choose to return 
% R_sw_ht_hi: upper bounds of the rates for the blocklengths specified in n_vals
% R_sw_ht_lo: lower bounds of the rates for the blocklengths specified in n_vals

    % Compute the optimal log_c coefficients from binary hypothesis testing, 
    % which serve as the start point for computing the parameters for 
    % composite hypothesis testing. 
    [~, log_c_binary] = p2p_ht_converse(a, n_vals, eps);

    % Compute log binomial coefficients
    log_b = binomial_coeff(max(n_vals));
    
    R_sw_ht_hi = zeros(n_vals,1);
    R_sw_ht_lo = zeros(n_vals,1);

    for ni = 1:length(n_vals)
        n = n_vals(ni);
        m = 0:n;
        log_Pr = m*log(a)+(n-m)*log(1/3*(1-a));  % log_Pr stores all the log joint probabilities in ascending order
        log_Pm = m*log(a+1/3*(1-a))+(n-m)*log(1/3*(1-a)+1/3*(1-a));  % log_Pm stores all the log marginal probabilities in ascending order

        % Initial points   
        x = log_c_binary(n) - 0.1;
        y = log_c_binary(n) - 0.0001;

        [beta1, beta2] = my_bisect(n,eps,x,log_Pr,log_Pm,log_b);
        bx = [beta1 beta2];
        diffx = beta1 - beta2 > 0;
        [beta1, beta2] = my_bisect(n,eps,y,log_Pr,log_Pm,log_b);
        by = [beta1 beta2];
        diffy = beta1 - beta2 > 0;

        % For catching special cases: boolean 1 -> zero exists, 0 -> no zero
        disp(strcat('n =',32,num2str(n),',',32,num2str(diffx ~= diffy)));  

        if diffx ~= diffy
            % Do bisection search if zero exists
            iter = 0;
            while abs(x-y) > 0.0001 && iter < 20 
                c = (x+y)/2;
                [beta1, beta2] = my_bisect(n,eps,c,log_Pr,log_Pm,log_b);
                if (beta1 - beta2 > 0) == diffx
                   x = c;     
                   bx = [beta1 beta2];
                else
                   y = c;
                   by = [beta1 beta2];
                end   
                iter = iter + 1;
            end
            [M,I] = min([max(bx) max(by)]);  % Minimum upper bound
            R_sw_ht_hi(ni) = M;
            lo = [min(bx) min(by)];
            R_sw_ht_lo(ni) = lo(I);          % A lower bound
        else
            % Scan a range of possible parameters for the best bound
            if (n < 100)
                c_val = log_c_binary(n)-0.1:0.0001:log_c_binary(n)-0.0001;
            else
                c_val = log_c_binary(n)-0.005:0.0001:log_c_binary(n)-0.0001;
            end
            beta1_val = zeros(length(c_val),1);
            beta2_val = zeros(length(c_val),1);
            for ci = 1:length(c_val)
                c = c_val(ci);
                [beta1, beta2] = my_bisect(n,eps,c,log_Pr,log_Pm,log_b);
                beta1_val(ci) = beta1;
                beta2_val(ci) = beta2;
            end
            R_sw_ht_hi(ni) = min(max(beta1_val,beta2_val));  % Minimum upper bound 
            R_sw_ht_lo(ni) = max(min(beta1_val,beta2_val));  % Maximum lower bound
        end
    end
    R = R_sw_ht_lo;
end



