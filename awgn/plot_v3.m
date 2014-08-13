function [Ns lb ub feinst gal wlfz] = plot_v3(A, epsil, hack, Ns)

if (nargin < 3) || (sum(size(hack)) == 0)
	hack = 0;
end

if(hack)
	disp('=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=');
	disp('=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=');
	disp('--> ATTENTION: test mode: kappa_inf() is used instead of kappa()  <--');
	disp('-->     However, estimated difference even for n = 10, A=1        <--');
	disp('-->     is 0.003 (for log M) or 3e-4 for rate.                    <--');
	disp('=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=');
	disp('=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=');
end

if (nargin < 4)  || (max(size(Ns)) == 0)
	% Ns = floor(10.^[1:.02:3.5]);
	Ns = 10:20:3010;
	%Ns = floor(linspace(1200, 1500, 10));
end

lb = [];
ub = [];
feinst = [];
gal = [];
wlfz = [];
taus = linspace(0,1,20).*epsil; taus = taus(3:end-2);

auto_hack = 0;
cycles_with_autohack = [];

for n = Ns
	%
	% TODO: need optimization over tau!!! THIS IS VERY IMPORTANT!!!
	%
	% For now I tested A = 1, epsil = 1e-3 :
	%  n=100: 	tau = epsil/2  -- almost optimal
	%  n=1000:	tau = epsil/5  -- the best, epsil/2 is .05% worse
	%
	disp(sprintf('--------------------- n=%d -----------------------', n));
	temp_lbs = []; 
	for tau = taus; 
		temp_lb = log2(kappa_inf(tau, A)) - betaq_up_v2(1-epsil+tau, n, A); 
		temp_lbs = [temp_lbs temp_lb]; 
	end;
	[bb ind] = max(temp_lbs);
	tau = taus(ind);
	if (hack) || (auto_hack)
		kap = log2(kappa_inf(tau, A));
	else
		kap = log2(kappa(tau, n, A));
	end
	clb = kap - betaq_up_v2(1 - epsil + tau, n, A);
	cub = -betaq_low_v2(1 - epsil, n, A);

	dif_bounds = 2*(clb - bb)/(clb + bb);
	
	disp(sprintf(	['     `-->  tally: \n'				...
			 '            best_tau = %.2f * epsil \n' 	...
			 '            dif ach bounds = %.2g %% \n'] ,	...
			 tau / epsil, 100*dif_bounds));
	disp(' ');

	if (abs(dif_bounds) < 1e-6) && (~hack) && (~auto_hack)
		auto_hack = 1;
		cycles_with_autohack = 0;
	end

	fst = feinstein_approx(n, epsil, A);

	feinst = [feinst fst/n];

	lb = [lb clb/n];
	ub = [ub cub/n];
	gal = [gal gallager_ach(n, epsil, A)/n];
	wlfz = [wlfz wolfowitz(n, epsil, A)/n];

	if (auto_hack)
		cycles_with_autohack = cycles_with_autohack + 1;
		if (cycles_with_autohack > 10)
			auto_hack = 0;
		end
	end
end

% generate 3 plots: comparing converses, comparing achievabilities, comparing approximation

Nsall = floor(linspace(Ns(1), Ns(end), 1000));
Ns_cap = [Ns(1) Ns(end)];
K = k_awgn(A, epsil);
Kapr = cap_awgn(A) + K ./ sqrt(Nsall);
Capr = cap_awgn(A) + 0 .* Ns_cap;

heps = -(1-epsil)*log2(1-epsil) - epsil*log2(epsil);
fano = (cap_awgn(A) + heps ./ Ns) ./ (1-epsil);

if A == 1
	ymax = 0.55;
elseif A == 10
	ymax = 3.5;
else
	ymax = 1.1*cap_awgn(A);
end;

%% Converses
figure;

plot(Ns_cap, Capr, 'r--', Ns, fano, 'k--', Ns, wlfz, 'k', Ns, ub, 'r');
xlabel('Blocklen, n'); ylabel('Rate, R'); ylim([0 ymax]);
title(sprintf('Converse bounds for AWGN (|X|^2=nP), SNR = %g dB, P_e = %g', 20*log10(A), epsil));
legend('Capacity', 'Fano', 'Wolfowitz', 'Converse', 'Location', 'Best');
grid on


%% Achievabilities
figure;
plot(Ns_cap, Capr, 'r--', Ns, lb, 'b', Ns, gal, 'b--', Ns, feinst, 'k--');
xlabel('Blocklen, n'); ylabel('Rate, R'); ylim([0 ymax]);
title(sprintf('Achievability bounds for AWGN, SNR = %g dB, P_e = %g', 20*log10(A), epsil));
legend('Capacity', 'New achievability', 'Gallager', 'Feinstein', 'Location', 'Best');
grid on

%% K-approximation
figure;
plot(Ns_cap, Capr, 'r--', Ns, ub, 'r', Nsall, Kapr, 'k-', Ns, lb, 'b', Ns, gal, 'b--');
xlabel('Blocklen, n'); ylabel('Rate, R'); ylim([0 ymax]);
title(sprintf('Bounds for AWGN, SNR = %g dB, P_e = %g', 20*log10(A), epsil));
legend('Capacity', 'Converse (|X|^2 = nP)',  'K-approximation', 'New achievability', ...
	'Gallager random coding', 'Location', 'Best' );
grid on


