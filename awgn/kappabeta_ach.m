function Lms = kappabeta_ach(Ns, epsil, P, hack)
%
% Compute achievability (lower) bound for log M. This is just a glue between kappa() and betaq_up_v2()
%

if (nargin < 4) || isempty(hack)
	hack = 1;
end


taus = linspace(0,1,40).*epsil; taus = taus(3:end-2);

Lms = [];

for n = Ns;
	disp(sprintf('kappabeta_ach(): n = %d', n));
	temp_lbs = []; 
	for tau = taus; 
		temp_lb = log2(kappa_inf(tau, P)) - betaq_up_v2(1-epsil+tau, n, P); 
		temp_lbs = [temp_lbs temp_lb]; 
	end;
	[bb ind] = max(temp_lbs);
	tau = taus(ind);
	if (hack) 
		kap = log2(kappa_inf(tau, P));
	else
		kap = log2(kappa(tau, n, P));
	end
	clb = kap - betaq_up_v2(1 - epsil + tau, n, P);

	Lms = [Lms clb];

end
