function f=normapx(Ns, delta, epsilon)
K=sqrt(delta*(1-delta))*log2((1-delta)/delta)*norminv(epsilon, 0, 1);
cap = 1 + delta.*log2(delta) + (1-delta).*log2(1-delta);
f = Ns.*cap + K.*sqrt(Ns) + log2(Ns)/2;
