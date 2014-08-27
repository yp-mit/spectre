function [C V] =  biawgn_stats(P)
% This function returns C and V for the BIAWGN channel.
% When using normal approximation don't forget the 1/2 log n factor.

% Experimental fact: Numerical integration fails for P > 19
if (P > 19)
	error('biawgn_stats: P too high!');
	return;
end

f = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log2(2) - log2(1+exp(-2*P+2*sqrt(P)*z)));

Zmin = -9; Zmax = 9;
if (abs(f(Zmin)) > 1e-10) || (abs(f(Zmax)) > 1e-10)
	disp(sprintf('ERROR: must increase range for P=%.4g', P));
	error('biawgn_stats: f-range');
end;

C = quadl(f, Zmin, Zmax);


g = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log2(2) - log2(1+exp(-2*P+2*sqrt(P)*z)) - C).^2;

if (abs(g(Zmin)) > 1e-10) || (abs(g(Zmax)) > 1e-10)
	disp(sprintf('ERROR: must increase range for P=%.4g', P));
	error('biawgn_stats: g-range');
end;

V = quadl(g, Zmin, Zmax);

