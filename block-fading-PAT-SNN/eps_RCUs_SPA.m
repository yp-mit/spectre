function eps_trial = eps_RCUs_SPA( i_s, n, L, r )

%allocate vectors for E0, E'0 and E''0
tau = linspace(0,1,500);
E0_tau=nan(size(tau));
E0_tau_1 = E0_tau;
E0_tau_2 = E0_tau;

for tt=1:length(tau) % Generation of derivatives of cumulant generating function
    % We substract the maximum in the exponent to avoid infinities
    E0_tau(tt)   =  -log(mean(exp(-tau(tt)*i_s)));
    E0_tau_1(tt)   = mean(exp(-tau(tt)*i_s).*i_s)/mean(exp(-tau(tt)*i_s));
    E0_tau_2(tt)   = (-mean(exp(-tau(tt)*i_s).*i_s.^2)*mean(exp(-tau(tt)*i_s))...
        +mean(exp(-tau(tt)*i_s).*i_s)^2)/mean(exp(-tau(tt)*i_s))^2;
end

V = -E0_tau_2;
A = real(sqrt((L*V)));

Pr_Q = exp(-L*(E0_tau-tau.*E0_tau_1)).*(exp(L/2*V.*tau.^2).*qfunc(tau.*A)...
    +exp((L/2*V.*(1-tau).^2)+log(qfunc((1-tau).*A))));

%Capturing possible issues:
Pr_Q(isnan(Pr_Q))=[]; %remove nans
Pr_Q(Pr_Q>1)=1; %if some probability larger than 1, replace with 1

[~, tau_indx] = max(E0_tau - tau*log(2^(r*n)-1)/L);
eps_trial = Pr_Q(tau_indx);


end