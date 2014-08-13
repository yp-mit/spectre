function lm = me_multin(Ns, epsil, A);
% Returns ldpc achevability data for a given (n,epsil,A)
cap = log2(1+A^2)/2;
Pdb = 20*log10(A);

lm = [];
for n = Ns;
	% these are absolute bounds on rate
	Rmin = 25/74; Rmax = 25/32;
	Rmax = min(Rmax, cap);
	% Step 1: find minimal acceptable rate
	Ptry = +Inf;
	Rnext = Rmin; Rstep = 0.1;
	while ~(Ptry <= Pdb)
		Rtry = Rnext;
		Ktry = floor(n*Rtry);
		[a b] = compute_ab(Ktry, Rtry);
		Ptry = (-norminv(epsil))/a + b;
		%% disp(sprintf('   Rtry = %.3g, Ktry = %d, Ptry = %.3g', Rtry, Ktry, Ptry));
		Rnext = Rnext + Rstep;
		if(Rnext > Rmax)
			if (Rstep == 0.1)
				Rnext = Rmin;
				Rstep = 0.01;
			else
				disp(sprintf(['me_multin(n = %d, epsil = %g, P = %.2g dB:' ...
						' could not find minrate'], n, epsil, Pdb));
				Rtry = Inf;
				break;
			end
		end
	end
	if(Rtry == Inf)
		lm = [lm NaN];
		continue;
	end
	Rmin = Rtry;
	Kmin = Ktry;
	Pmin = Ptry;
	% Step 2: find maximal acceptable rate

	Ptry = -Inf;
	Rnext = Rmax; Rstep = 0.1;
	while ~(Ptry >= Pdb)
		Rtry = Rnext;
		Ktry = ceil(n*Rtry);
		[a b] = compute_ab(Ktry, Rtry);
		Ptry = (-norminv(epsil))/a + b;
		Rnext = Rnext - Rstep;
		if(Rnext < Rmin)
			if (Rstep == 0.1)
				Rnext = Rmax;
				Rstep = 0.01;
			else
				disp(sprintf(['me_multin(n = %d, epsil = %g, P = %.2g dB:' ...
						' could not find maxrate'], n, epsil, Pdb));
				Rtry = Inf;
				break;
			end
		end
	end
	if(Rtry == Inf)
		lm = [lm NaN];
		continue;
	end
	Rmax = Rtry;
	Kmax = Ktry;
	Pmax = Ptry;
	% Step 3: proceed by binary division

	disp(sprintf(['me_multin(n = %d, eps = %g, P = %.2g dB: R = (%.3g, %.3g), '...
			'K = (%d, %d), P = (%.3g, %.3g)'], ...
			n, epsil, Pdb, Rmin, Rmax, Kmin, Kmax, Pmin, Pmax));

	while (Kmax - Kmin) > 1
		Rtry = (Rmin + Rmax)/2;
		Ktry = floor(Rtry*n);
		[a b] = compute_ab(Ktry, Rtry);
		Ptry = (-norminv(epsil))/a + b;
		%% disp(sprintf('   Rtry = %.3g, Ktry = %d, Ptry = %.3g', Rtry, Ktry, Ptry));
		if(Ptry == NaN) break; end
		if(Ptry > Pdb)
			Rmax = Rtry; Kmax = Ktry;
		else
			Rmin = Rtry; Kmin = Ktry;
		end
	end
	reason = 'PREC OK:';
	if(Ptry == NaN)
		reason = 'NAN BREAK:';
	end
	disp(sprintf('     %s lm = %d', reason, Kmin));
	lm = [lm Kmin];
end
