function Lms = normapx_awgn(Ns, epsil, P);

C = cap_awgn(P);
K = sqrt(P)./sqrt(2) .* sqrt(P + 2)./(P + 1) *log2(exp(1)) * norminv(epsil, 0, 1);
Lms = C * Ns + K * sqrt(Ns) + log2(Ns)/2;
