function lp = ncx2log(qq0,n,delta)

% returns log(ncx2cdf(qq0,n,delta)) but is tailored for v = n, delta ~ C*n, x = K n, 
% n->infty
%
% results are consistent for lp > -1240. Below that everything blows up fast


accum = -Inf;

ii = 0;

flag_peak = 0;

hdelta = delta/2;

while 1
	% for ii >= 25 precision of Gosper's formula for n! is < 1e-5
	% and is below 1e-8 for ii's of interest
	if ii < 25
		pp = poisspdf(ii, hdelta);
		if pp == 0
			lpois = -Inf;
		else
			lpois = log(pp);
		end
	else
		
		lpois = -hdelta + log(hdelta) * ii - ... % ii !
			( .5 * log((2*ii+1/3)*pi) + ii * log(ii) - ii );
	end

	if(mod(ii, 10) == 0)
		GG = gammainc(qq0/2, n/2 + ii);
	else
		incr = exp(-qq0/2 + log(qq0/2)*(n/2+ii-1) - gammaln(n/2 + ii));
		if ( abs(GG - incr) < 1e-6 * GG ) || (GG < incr)
			% instability! Recalculate
			GG = gammainc(qq0/2, n/2 + ii);
		else
			GG = GG - incr;
		end
	end
	
	if GG == 0
		lGG = -Inf;
	else
		lGG = log(GG);
	end	
	lt = lpois + lGG;

	if(lt > -Inf)
		
		flag_peak = 1;

		if(lt > accum)
			accum = lt + log(1+exp(accum-lt));
		else
			accum = accum + log(1+exp(lt - accum));
		end
	end

	ii = ii + 1;

	if ( ii > 5*delta ) && (flag_peak == 0)
		%disp('!!! !!! Warning !!! !!!: ncx2log(): returing -Inf');
		break;
	end

	if (ii > 5*delta) && (flag_peak == 1)
		disp('Warning: ncx2log(): something strange going on... (might hang)');
	end

	if (lt == -Inf) && (flag_peak == 1)
		% Ok, we passed the peak and lt = -Inf, terminate!
		break;
	end

end

lp = accum;
