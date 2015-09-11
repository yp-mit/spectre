function logM = MIMO_csir_normapx(n,epsilon,n_t,n_r,T,P,real_or_complex)

[C,V] = capacity_and_dispersion(n_t,n_r,T,P,real_or_complex);
logM = n*T*C - sqrt(n*T*V)*qfuncinv(epsilon);