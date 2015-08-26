function En0 = energy_awgn_ach(Lms, epsil);
% I use the following bound on the orig. ach-ty bound:
% denote by x* anything satisfying M Q(x*) > 1
%
% then
% epsil <= \Phi(x*-sqrt(E)) + M*exp(-x*^2/2 - (x*-sqrt(E))^2/2)/(1+x*^2)/2pi
%
% Warning: fzero fails for lm > 1e4

En0 = [];
for lm = Lms; 
%	disp(sprintf('Starting lm=%g\n', lm));
	% to guess x_* use the lower bound on Q-function:
	% Q(x) > x/(1+x^2) e^(-x^2/2) / sqrt(2pi)
	x0 = sqrt(2*lm*log(2));
	if (x0 < 11)
		func1 = @(x) normcdf(-x)*2^lm - 1;
	else
		func1 = @(x) x./(1+x.^2) .*exp(-x.^2/2 + lm*log(2))/sqrt(2*pi) - 1;
	end
	xs = fzero(func1, x0);
	%keyboard;
	func2 = @(E) normcdf(xs-sqrt(E)) + exp(-xs.^2/2 - (xs-sqrt(E)).^2/2+lm*log(2)-log(1+xs.^2))/(2*pi);

	E0 = (xs-norminv(epsil))^2;
	% Es = linspace(E0*.8, E0*1.2); figure(1); plot(Es, func2(Es)); pause;

	E = fzero(@(E) func2(E) - epsil, E0);
	En0 = [En0 E/2];
end
