function P = optpower(epsil, n, R);
% compute the smallest required power to communicate at rate R with p.e. epsil and blocklength n
K = @(A) A./sqrt(2) .* sqrt(A.^2 + 2)./(A.^2 + 1) *log2(exp(1)) * norminv(epsil, 0, 1);

normapx = @(A) 1/2 * log2(1 + A.^2) + K(A)/sqrt(n) + 1/2 *log2(n)/n;
A0 = sqrt(2^(2*R) - 1);
A1 = fzero(@(A) normapx(A) - R, A0);
P = A1^2;
