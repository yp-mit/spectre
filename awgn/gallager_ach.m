function lm = gallager_ach(n, epsil, P)
% returns log M achievable by gallager's random coding.
%
% Note: In the case of uncostrained input gallager's bound can be written as
%
%  P_e <= M^s   V(s)^n
%
% Where V(s) = (   Sum_y [ Sum_x Q(x) W(y|x)^{1/1+s} ]^(1+s)   )
%
% Thus we could rewrite the bound as a lower bound on $M$ directly.
%
% For AWGN we can do the same in the case when we do not restrict codewords to a narrow shell 
% around the power sphere. I did that and here is the result:
%
%   gallager_ach(3000, 1e-6, 1) = 1225
%   -betaq_low_v2(1-1e-6, 3000, 1) = 1287
%  
%   new_gallager() = 1206 --- much worse than good gallager!
%
% Here's how to compute new_gallager:
% mu = chi2cdf(n,n); 
% f = @(x) -( (1+x)./x * log2(mu) + n/2 .* log2( 1 + A^2 ./ (1+x)) + log2(epsil/2) ./ x )
% [x_best f_best] = fminbnd(f, 0, 1); new_gal = -f - 1;
%
% Note epsil/2 and also that we should subtract one bit.

R_up = cap_awgn(P);
R_down = 0;

% Check if for this n gallager works at all
Pe = gallager_pe(P, R_down, n);
if Pe > epsil;
    lm = 0;
    return;
end

precision = 1e-3 * R_up / n;

while 1;
	if ((R_up - R_down) < precision)
		break;
	end
	R = (R_up + R_down) / 2;

	Pe = gallager_pe(P, R, n);
	if( Pe < epsil)
		R_down = R;
	else
		R_up = R;
	end
end

lm = n * R_down;


function [pe Er] = gallager_pe(P, R, n, deltap_force)
%
% This function returns upper bound (achievability) on P_e using
% Gallager's random coding error exponents
%
% TODO: use expurgated exponent for low rates!
%

% conversion A->P. Old versions are all in terms of ``amplitude'' A.
A = sqrt(P);

if (nargin < 4)
	deltap_force = [];
end

if (R == 0)
	pe = 0;
	return;
end

Rcr = 1/2 * log2( 1/2 + A^2/4 + 1/2 * sqrt(1+A^4/4) );

% convert rate to nats
R = R / log2(exp(1));
Rcrn = Rcr / log2(exp(1));


if (R >= Rcrn)
	beta = exp(2*R);
	ro = A^2/(2*beta) * (1 + sqrt(1+ 4*beta/A^2/(beta-1) ) ) - 1;

	Er = A^2/(4*beta) * ...
		( (beta+1) - (beta-1)*sqrt(1 + 4*beta/A^2/(beta-1)) ) + ...
		1/2 * log(  beta - A^2*(beta-1)/2 * ( sqrt(1 + 4*beta/A^2/(beta-1)) - 1 ) );
else
	beta = 1/2 * (1 + A^2/2 + sqrt(1+A^4/4) );
	ro = 1;

	Er = 1 - beta + A^2/2 + 1/2 * log(beta - A^2/2) + 1/2*log(beta) - R;
end;


s = ro*A^2/(2 * (1+ro)^2 * beta);

%%
%% Now instead of maximizing over delta's we choose delta = 1/s  (see Gallager 7.4.39)
%% This is suboptimal, but let it be...
%% 

if (max(size(deltap_force)) > 0)
	deltap = deltap_force * n;
else
	if (s == 0) 
		deltap = n; 
	else
		deltap = min(1/s, n);
	end;
end

%disp(sprintf('ro = %g, deltap = %g, s = %g', ro, deltap, s));

ccdfs = chi2cdf([n-deltap n], n);

mu = ccdfs(2) - ccdfs(1);
if (mu < 1e-10 * ccdfs(1))
	mu_new = quad(@(x)( chi2pdf(x, n) ), n - deltap, n);
	disp(sprintf( [	'--- gallager_pe: computing mu using quad()\n' ...
			'    `--> mu_new = %.3g, mu_old = %.3g (diff = %.2f%%)'], ...
			mu_new, mu, 200*abs(mu_new - mu)/(mu_new + mu)));
	mu = mu_new;
end;

multip = 2*exp(s*deltap)/mu;
multip_approx = ro*A^2*exp(1)*sqrt(4*pi*n)/(1+ro)^2 /beta;

%	disp(sprintf( [	'--- gallager_pe: multiplier comparsion\n' ...
%			'    `--> multip = %.3g, approx = %.3g (diff = %.2f%%)'], ...
%			multip, multip_approx, 100*abs(multip - multip_approx)/(multip)));


pe = multip*exp(-n*Er);

if(pe > 1) pe = 1; end;

