function log_b = binomial_coeff(n_max)
    % Compute log binomial coefficients
    log_b = {};
    for n = 0:n_max
        Ks = 0:n;
        log_b{n+1} = gammaln(n+1) - gammaln(Ks+1) - gammaln(n-Ks+1);
    end
end