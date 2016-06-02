function [Lms Lms_avg] = dt_ach(Ns, delta, epsil);
% Compute DT achievability bound:
%
%   P_e <= mdist(N, delta, log(M/2)).
%
% The same upper-bound holds for MAXIMAL prob. of error provided that M=2^k (random linear code).
%
% NOTE: P_und is the average undetected prob. of error %

shutup = 0;

if max(size(Ns,2)) > 1
        shutup = 1;
else
	disp(sprintf(' -- dt_ach(n=%d, delta = %g, epsil = %g): ', Ns(1), delta, epsil));
end

Lms = zeros(size(Ns));
Lms_avg = zeros(size(Ns));

for idx = 1:size(Ns,2)
	n = Ns(idx);

	% Our bound is useless?
	if(mdist(n, delta, 0) > log2(epsil))
		Lms(idx) = 0; Lms_avg(idx) = 0;
		continue;
	end

	gam0 = 0;	% Can always choose 1 codeword)
	gam1 = n;	% Can not choose more than 2^n

	% use normal approximation for the first shot
	gam_cur = normapx(n, delta, epsil);
	step = 1;

	% Repeat until our gam's do not surround any integer
	while( floor(gam0) < floor(gam1) )
		[test pund] = mdist(n, delta, gam_cur);
		if(test < log2(epsil))
			if ~shutup
				disp(sprintf(['   Step %d. gam = [%.1f, %.1f], tested: %.1f, test=%.3f, '...
						'Heading UP'],...
						step, gam0, gam1, gam_cur, 2^test));
			end
			gam0 = gam_cur;
		else
			if ~shutup
				disp(sprintf(['   Step %d. gam = [%.1f, %.1f], tested: %.1f, test=%.3f, '...
						'Heading DOWN'],...
						step, gam0, gam1, gam_cur, 2^test));
			end
			gam1 = gam_cur;
		end
		step = step + 1;
		gam_cur = (gam0 + gam1)/2;
	end
	Lms(idx) = floor(gam0) + 1;
	Lms_avg(idx) = gam0;
	if ~shutup
		disp(sprintf('      ====> log M >= %g', lm(idx)));
		disp(sprintf('            P_und <= %g (on avg)', 2^pund));
	end
end


function [f g] = mdist(n, delta, Gammas)
% returns the f=log2 of value of P[i(X,Y) <= \gamma] + 2^\gamma P[i(X,\bar Y) > \gamma] for BSC
% and g = log2 of the product

A = 1 + log2(1-delta); B = log2((1-delta)/delta);
f = zeros(1,size(Gammas,2));
g = f;
Ks = 0:n;

% This method computes
%
%    E 2^{i(X, \bar Y)} \wedge 2^gamma
%
terms0 = (gammaln(n+1) - gammaln(Ks+1) - gammaln(n-Ks+1))/log(2);

infdens = n*A-B*Ks;

for idx = 1:size(Gammas,2);
	gamma = Gammas(idx);


	terms = terms0 + min(infdens, gamma)-n;

	f(idx) = sumlog2(terms);

	terms1 = terms(infdens<=gamma);
	terms2 = terms(infdens>gamma);

	g(idx) = sumlog2(terms1)+sumlog2(terms2);
end
