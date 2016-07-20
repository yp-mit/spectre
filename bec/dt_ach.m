function [lm lm_avg] = achiev2(Ns, delta, epsil);
% Compute DT achievability bound. lm is the bound for maximal Pe. lm_avg is the bound for average Pe.
%
% For BEC achievability bound is (Pe-max!)
%
%   epsil <= mdist(N, delta, log(M)).
%
% In fact, using random-linear-code trick we can do better (again Pe-max!)
%   
%   epsil <= mdist(N, delta, ceil(log(M/2)))
%
%

shutup = 0;

if max(size(Ns,2)) > 1
        shutup = 1;
else
	disp(sprintf(' -- achiev(n=%d, delta = %g, epsil = %g): ', Ns(1), delta, epsil));
end

lm = zeros(size(Ns)); lm_avg = lm;

for idx = 1:size(Ns,2)
	n = Ns(idx);

	% Our bound is useless?
	if(mdist(n, delta, 0) > log2(epsil))
		lm(idx) = 0;
		continue;
	end

	gam0 = 0;	% Can always choose 1 codeword)
	gam1 = n;	% Can not choose more than 2^n

	% use normal approximation for the first shot
	gam_cur = normapx(n, delta, epsil);
	step = 1;

	% Repeat until our gam's are not close 
	while( (gam1-gam0) > 1e-4 )
		test = mdist(n, delta, gam_cur);
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
	lm_avg(idx) = gam0+1;
	lm(idx) = floor(gam0+1);
	if ~shutup
		disp(sprintf('      ====> log M >= %g', lm(idx)));
	end
end
