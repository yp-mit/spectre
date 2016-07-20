function Lms=normapx(Ns, delta, epsilon)
K=sqrt(delta*(1-delta))*norminv(epsilon, 0, 1);
Lms = Ns.*(1-delta) + K.*sqrt(Ns);
