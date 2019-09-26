function log_x = my_logdiffexp(a,b)
% Self written function to compute the log difference of the exponentials of a, b.
% This function avoids overflow in computing the exponentials.
% Input: 
% a, b: two values
% Output:
% log_x: log of the difference of the exponentials of a, b
    if a > b
        log_x = a + log(1-exp(b-a));
    else
        log_x = -inf;
    end
end