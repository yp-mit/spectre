function log_y = my_logsumexp(vec)
% Self written function to compute the log sum of the exponentials of each element in the input vector.
% This function avoids overflow in computing the exponentials.
% Input: 
% vec: a vector
% Output:
% log_y: log of the sum of the exponentials
    s = max(vec);
    if s ~= -inf
        log_y = s + log(sum(exp(vec-s)));
    else
        log_y = -inf;
    end
end