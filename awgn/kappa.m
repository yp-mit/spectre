function f = kappa(tau, n, P)
% computes:
%	kappa(tau) = min Q(B), s.t. P(B) \ge tau 
%
%	P, Q defined in dpdq.m
%
%
% TODO: kappa(1e-7, 93, 10) doesn't work...  Can it be fixed?

% conversion A->P. Old versions are all in terms of ``amplitude'' A.
A = sqrt(P);

posval = @(x) (abs(x) + x)./2;
negval = @(x) (abs(x) - x)./2;

if(tau <= 0 || tau > 1)
	disp(sprintf('ERROR: kappa(): tau = %g, bad boy!', tau));
	error('kappa')
end
midpoint = (1+A^2)*n;
xs = linspace(.25*midpoint, 2.5*midpoint, 10000);
ds = dpdq(xs, n, A);

[y in] = max(ds);

%gammas = linspace(ds(1), y, 1000);

gamma_b = ds(1);
gamma_u = y;

iter = 1;

%
% Iterate over trial gammas
% 
while 1;
	
	gamma = (gamma_b + gamma_u)/2;
	%
	% Step 1. Compute x_1, x_2 
	%
	inx1 = find(posval(ds(1:in) - gamma),1);
	if (inx1 == 1) 
		disp(sprintf(['ERROR: kappa(tau = %g, n = %d, A = %g): inx1 = %d, in = %d, iter = %d --' ...
			'can not proceed.'], tau, n, A, inx1, in, iter));
		error('kappa', 'k');
	end

	xa = xs(inx1-1);
	xb = xs(inx1);

	x1 = xa + (gamma - ds(inx1-1))/(ds(inx1) - ds(inx1-1)) * (xb - xa);

	inx2 = find(negval(ds(in:end) - gamma), 1);
	if(sum(size(inx2)) == 0)  || (inx2 == 1)
		disp(sprintf(['ERROR: kappa(tau = %g, n = %d, A = %g): inx2 = %d, in = %d --' ...
			'can not proceed.'], tau, n, A, inx2, in));
		error('kappa', 'k2');
	end

	inx2 = in + inx2 - 1;

	xa = xs(inx2 - 1);
	xb = xs(inx2);

	x2 = xa + ( ds(inx2 - 1) - gamma ) / ( ds(inx2 - 1) - ds(inx2) ) * (xb - xa);

	%
	% Step 2. Compute resulting tau: two methods pdf integration and 
	%

	P = ncx2cdf([x1 x2], n, n*A^2);
	curtau1 = P(2) - P(1);
	
	curtau2 = quad(@(x)ncx2pdf(x, n, n*A^2), x1, x2);

	if (2*(curtau2 - curtau1)/(curtau2+curtau1) > 1e-5)
		disp( sprintf( ...
			'WARNING: kappa(n=%d, A=%g): gamma = %g, delta (ctau1, ctau2) = %g %%', ...
			n, A, gamma, 200*(curtau2 - curtau1)/(curtau2+curtau1)));
	end

	curtau = (curtau2 + curtau1)/2;

	if (abs(tau - curtau) / tau) < 1e-6
		break;
	end

	if (curtau > tau) 
		gamma_b = gamma;
	else
		gamma_u = gamma;
	end

	iter = iter + 1;
	
end

%
% Again, compute kappa using two methods
%
	
P = gamcdf([x1 x2], n/2, 2*(1+A^2));
kappa1 = P(2) - P(1);

kappa2 = quad(@(x)gampdf(x, n/2, 2*(1+A^2)), x1, x2);

if (2*(kappa2 - kappa1)/(kappa2+kappa1) > 1e-5)
	disp( sprintf( ...
		'WARNING: kappa(n=%d, A=%g): gamma = %g, delta (kappa1, kappa2) = %g %%', ...
		n, A, gamma, 200*(kappa2 - kappa1)/(kappa2+kappa1)));
end
kappa = (kappa2 + kappa1)/2;

kapinf = kappa_inf(tau, P);

delta = 200 * (kappa - kapinf) / (kappa + kapinf);
	
disp(sprintf( ...
	['kappa(tau = %g, n = %d, A = %g): ' 			...
   		'iter = %d, gamma = %g, x_1 = %g, x_2 = %g\n' 	...
	 '        kappa = %g, kapinf = %g, delta = %g %%'], 	...
	tau, n, A, iter, gamma, x1, x2, kappa, kapinf, delta ));

f = kappa;

function f = dpdq(x, n, A)
%	computes dP/dQ derivative at x for
%		P ~ sum (Z_i + A)^2
%		Q ~ sum (1+A^2) Z_i^2
%
%	This is crucial part for computing kappa(tau)

P = ncx2pdf(x, n, n*A^2);
Q = gampdf(x, n/2, 2*(1+A^2));

% I know that if Q is very small then P is even smaller
Q = Q + (10^10 .* (Q == 0));

f = P ./ Q;

%if Q < 100*eps
%	if P < 100*eps
%		disp(sprintf('!!! ERROR: dpdq(x = %g, n = %d, A = %g): Q = %g, P= %g ', ...
%			x, n, A, Q, P));
%		f = P/Q;
%	else
%		f = +Inf;
%	end
%else
%	f = P / Q;
%end

