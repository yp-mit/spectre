function K = k_awgn(A, epsil)
K = A./sqrt(2) .* sqrt(A.^2 + 2)./(A.^2 + 1) *log2(exp(1)) * norminv(epsil, 0, 1);
