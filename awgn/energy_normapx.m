function lm = normapx_nofb(en0, epsil)
% Compute normapx on log m^*(E/N_0, \epsilon)

loge = log2(exp(1));
lm = en0 * loge + loge * sqrt(2*en0) * norminv(epsil) + log2(en0)/2;
%lm = en0 * loge + loge * sqrt(2*en0) * norminv(epsil) + log2(2*en0)/2;
