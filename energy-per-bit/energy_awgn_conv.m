function lm = energy_awgn_conv(en0, epsil)
% Compute upper bound on log m^*(E/N_0, \epsilon)

% Use the lower bound on Q-function:
% Q(x) > x/(1+x^2) e^(-x^2/2) / sqrt(2pi)
argu = sqrt(2*en0) + norminv(epsil);
lm = -(log2(argu./(1+argu.^2)) - log2(2*pi)/2 - argu.^2 * log2(exp(1))/2);
lm(argu < 11) = -log2(normcdf(-argu(argu<11)));
