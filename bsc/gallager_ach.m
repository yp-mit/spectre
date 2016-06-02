function Lms = gallager_ach(Ns, delta, epsil)
% computes gallager's achievability bound for BSC 
%
% We use it in the form 
%
%     P_e \le M^s  (   Sum_y [ Sum_x Q(x) W(y|x)^{1/1+s} ]^(1+s)   )^n
%
% Note: we then conclude that from (M, P_e) code with AVG P_e 
% we can get (M/2, 2 P_e) code with MAXIMAL. Thus we must take P_e = epsil/2 and then 
% subtract 1 bit from gallager's result
% 

Lms = [];

for n = Ns;

	f = @(x) -( log2(epsil/2) ./ x + n - n .* (1+x) ./ x .* log2( delta.^(1./(1+x)) + (1-delta).^(1./(1+x))));

	[lambda f] = fminbnd(f, 0, 1);

	% Cut half of the codewords

	logm = -f-1;

	disp(sprintf('--- Gallager: best_lambda = %g, log M >= %g', lambda, logm));

	if ( lambda <= 0 ) || (lambda >= 1)
		disp(sprintf('------- WARNING: lambda on the boundary! params: n=%d, delta=%g, epsil=%g', n, delta, epsil));
	end
	Lms = [Lms logm];
end
