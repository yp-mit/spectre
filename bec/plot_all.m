function [Ns ach_dt ach_gal conv normapx_val] = plot_all(delta, epsil);

if (nargin < 1) || isempty(delta)
	delta = 0.5;
end
if (nargin < 2) || isempty(epsil)
	epsil = 1e-3;
end

Ns = 20:20:1000;
Ns_cap = [Ns(1) Ns(end)];
cap = 1 - delta;
Capr = cap + Ns_cap*0;


normapx_val = normapx(Ns, delta, epsil);
conv = converse(Ns, delta,epsil);
ach_gal = gallager_ach(Ns, delta, epsil);
ach_dt = dt_ach(Ns, delta, epsil);


%% normal approximation
figure;
plot(Ns_cap, Capr, 'r--', Ns, conv./Ns, 'r', Ns, normapx_val./Ns, 'k-',  Ns, ach_dt./Ns, 'b', ...
	Ns, ach_gal./Ns, 'b-.');
xlabel('Blocklen, n'); ylabel('Rate, R'); ylim([0 1.05*cap]);
title(sprintf('Bounds for the BEC(%g), P_{e,max} = %g', delta, epsil));
legend('Capacity', 'Converse',  'Normal approximation', 'DT achievability', 'Gallager achievability', 'Location', 'SouthEast');
grid on


