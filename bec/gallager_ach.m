function Lms = gallager_ach(Ns, delta, epsil)
% computes gallager's achievability bound for BEC
%
% We use it in the form 
%
%     P_e \le M^s  (   Sum_y [ Sum_x Q(x) W(y|x)^{1/1+s} ]^(1+s)   )^n
%
%
% Note: we then conclude that from (M, P_e) code with AVG P_e 
% we can get (M/2, 2 P_e) code with MAXIMAL. Thus we must take P_e = epsil/2 and then 
% subtract 1 bit from gallager's result
%
% For BEC E_0(s) = -log2(delta + 2^-s (1-delta));

Lms = [];

for n = Ns;

	f = @(x) -( log2(epsil/2) -n*log2(delta + 2.^(-x) * (1-delta)) )./x;

	[lambda f] = fminbnd(f, 0, 1);

	% We drop one bit as discussed above.
	logm = -f-1;

	disp(sprintf('--- Gallager: n = %d, best_lambda = %g, log M >= %g', n, lambda, logm));

	if ( lambda <= 0 ) || (lambda >= 1)
		disp('------- WARNING: lambda on the boundary!');
	end

	Lms = [Lms logm];
end
