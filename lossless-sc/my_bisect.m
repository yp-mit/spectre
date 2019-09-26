function [beta1, beta2] = my_bisect(n, eps, c, log_Pr, log_Pm, log_b)
% This is a helper function called by sw_ht_converse in computing the
% parameters for the Slepian-Wolf hypothesis testing converse.

% Input: 
% n: blocklength
% eps: target error
% c: parameter in hypothesis testing
% log_Pr: joint probabilities on log scale
% log_Pm: marginal probabilities on log scale
% log_b: binomial coefficients on log scale

    log_c = n*c;
    % Screen all type classes
    log_A = [];
    for m = 0:n
        for j = 0:n-m
            for k = 0:n-m-j
                log_A_tempt = my_logdiffexp(log_Pr(m+1),log_c) - my_logsumexp([log_Pm(m+j+1) log_Pm(m+k+1)]);
                log_A = [log_A; [log_A_tempt m j k]];   
            end
        end
    end
    log_A = sortrows(log_A,'descend');
    log_A(log_A(:,1) == -inf,:) = [];

    % Add up the probabilities of type classes in descending order until the total exceeds 1 - eps
    i = 1;
    m = log_A(i,2); j = log_A(i,3); k = log_A(i,4);
    log_sum_P = log_Pr(m+1) + log_b{n+1}(m+1)+log_b{n-m+1}(j+1)+log_b{n-m-j+1}(k+1);
    log_U = [];
    log_UP = [];
    while exp(log_sum_P) < 1 - eps && i < size(log_A,1)
        i = i + 1;
        m = log_A(i,2); j = log_A(i,3); k = log_A(i,4);
        log_U_tempt = log_b{n+1}(m+1)+log_b{n-m+1}(j+1)+log_b{n-m-j+1}(k+1);
        log_P_tempt = log_Pr(m+1) + log_U_tempt; 
        log_UP_tempt = log_Pm(m+j+1) + log_U_tempt;

        log_sum_P = my_logsumexp([log_sum_P log_P_tempt]); 
        log_U = [log_U; log_U_tempt];
        log_UP = [log_UP; log_UP_tempt];
    end

    if exp(log_sum_P) >= 1 - eps
        thresh_index = find(log_A(:,1) == log_A(i,1))';
        log_B_P_cut = [];
        log_B_P = [];
        log_B_U = [];
        log_B_UP = [];
        for t = thresh_index
            m = log_A(t,2); j = log_A(t,3); k = log_A(t,4);
            log_B_U_tempt = log_b{n+1}(m+1)+log_b{n-m+1}(j+1)+log_b{n-m-j+1}(k+1);
            if t <= i
                log_B_P_cut = [log_B_P_cut; log_Pr(m+1) + log_B_U_tempt]; 
            end
            log_B_P = [log_B_P; log_Pr(m+1) + log_B_U_tempt];  
            log_B_U = [log_B_U; log_B_U_tempt];
            log_B_UP = [log_B_UP; log_Pm(m+j+1) + log_B_U_tempt];
        end    
        gap = 1 - eps - (exp(log_sum_P) - exp(my_logsumexp(log_B_P_cut))); 
        log_lambda = log(gap) - my_logsumexp(log_B_P);  

        % Compute log_beta_2
        log_U = log_U(1:thresh_index(1)-1);
        U1 = my_logsumexp(log_U);   
        U2 = log_lambda + my_logsumexp(log_B_U);
        log_beta_2 = my_logsumexp([U1 U2]);

        % Compute log_beta_1
        log_UP = log_UP(1:thresh_index(1)-1);
        UP1 = my_logsumexp(log_UP);   
        UP2 = log_lambda + my_logsumexp(log_B_UP);
        log_beta_1 = my_logsumexp([UP1 UP2]);

        beta1 = 1/n*2*log_beta_1;
        beta2 = 1/n*log_beta_2;
    else
        % Special values to indicate error.
        beta1 = 100;
        beta2 = 0;
    end
end