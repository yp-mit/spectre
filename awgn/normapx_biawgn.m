function Lms = normapx_biawgn(Ns, epsil, P);
[C V] = biawgn_stats(P);
K = sqrt(V) * norminv(epsil, 0, 1);
Lms = C * Ns + K * sqrt(Ns) + log2(Ns)/2;
