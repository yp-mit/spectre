function [lm theta1] = shannon_ach2(n, epsil, P)
%	returns Shannon's achievability bound for fixed-power code of SNR A^2
%	(maximal prob. of error)
%
%	The difference from shannon_ach() is that we optimize ``purging'' by maximizing \tau:
%	(1-\tau)M_S(n, \tau * \epsilon)
%	
%	This is done approximately with the use of normal approximation.
%
% Note: qn_low, qn_up and nctcdf0 are helper functions defined below.

% conversion A->P. Old versions are all in terms of ``amplitude'' A.
A = sqrt(P);

V = A./sqrt(2) .* sqrt(A.^2 + 2)./(A.^2 + 1) *log2(exp(1));

func = @(x) -log2(1-x) - sqrt(n*V)* norminv(x*epsil);
tau = fminbnd(func, 0, 1-1e-6);


% Step 1. Find minimum theta1 such that 
% 	qn_up(theta1) <= epsil ;

% we must find some high theta for which qn_up() <= epsil
% atan(sqrt(n)) is the maximum possible one, so this one is pretty safe too.

theta_hi = atan(sqrt(n-1));

if (qn_up(n, epsil, A, theta_hi) > tau*epsil)
	 disp('--- shannon_ach(): can not find suitable theta, Log M >= 0');
	 theta1 = pi/2;
	 lm = 0;
	 return;
end

theta_lo = 0;

% This should give log M with relative prec. 10^-3 (approx)
logm_precision = 1e-3/2;
iter = 0;

fn_hi = +Inf;
fn_lo = -Inf;

F_hi = 0;
F_lo = 1;

while (true);
	theta = (theta_hi + theta_lo)/2;

	F = qn_up(n, epsil, A, theta);

	if (F > epsil)
		theta_lo = theta;
		F_lo = F;
	else
		theta_hi = theta;
		F_hi = F;
	end

	iter = iter + 1;
	
%	disp(sprintf('-- iter = %d, theta = [%.5g; %.5g]; (F-eps/2) = %.5g * eps', ...
%			iter, theta_lo, theta_hi, (F-tau*epsil)/epsil));

	if(iter > 200)
		error('ERROR: binary search does not converge...');
		return;
	end;
	
	if (theta_lo > 0) && (theta_hi < pi/2)
		fn_lo = logfntheta(n, theta_lo);
		fn_hi = logfntheta(n, theta_hi);

		%disp(sprintf('      fn = [%.5g, %.5g], deltafn = %.5f (%.4f %%)', ...
		%	fn_lo, fn_hi, fn_hi - fn_lo, 200*(fn_hi - fn_lo)/abs(fn_hi + fn_lo) ));

		if( (fn_hi - fn_lo) < logm_precision * abs(fn_lo+fn_hi) )
			break;
		end
	end

end

theta1 = theta_hi;

% Step 2. Calculate lower bound on log M using fn(theta)
%
% Log M \ge -log(fn(theta)));
% fn(theta) = gamma(n/2 + 1)/gamma(n/2 + 1/2) / n / sqrt(pi) * sin(theta)^(n-1)/cos(theta)
%

%sin_th = sin(theta1);
%tg_th = tan(theta1);

fn_th = fn_hi;

lm = -fn_th + log2(1-tau) ;


disp(sprintf('--- shannon_ach(n=%d, eps=%.3g, P=%g): log M >= %.4g   (%d iters)', ...
			n, epsil, P, lm, iter));
disp(sprintf('         precision: tau*eps - Qup(theta) = %g ', ...
			tau*epsil - F_hi));
disp(sprintf('         logM precision: %.4f (%.4f %%)', ...
			fn_hi - fn_lo, 100*(fn_hi - fn_lo) / lm));
disp(sprintf('         theta_opt = %.8g, tau_opt = %.6g', theta1, tau));



function q = qn_low(n, epsil, A, theta)
% This function is a Shannon's lower bound on P_e,opt of fixed-power code
%
% We need epsil only for compat with qn_up()
%

 q = 0*theta;

 q( theta<= 0 ) = 1;
 q( theta >= pi/2) = 0;

 idx = (theta > 0) & (theta < pi/2);
 
 nu = n-1;
 delta = sqrt(n) * A;
 q(idx) = nctcdf0(sqrt(n-1)*cot(theta(idx)), nu, delta);

function f = qn_up(n, epsil, A, theta1)
%	This function is an upper-bound on the P_e,opt of best fixed-power code of
%	\approx 1/f_n(theta1) codewords.

if (theta1 == 0)
	f = 1;
	return;
end

multip = cot(theta1)/(1-tan(theta1)^2/n) * (n-1);

if(multip < 0)
	error('ERROR: sorry, for theta1 so large, this bound does not work');
	return;
end


f_int = @(theta) qn_low(n, epsil, A, theta) .* (sin(theta) ./ sin(theta1)).^(n-2);

integral = quadl(f_int, 0, theta1, 1e-6*epsil);

f = multip * integral;


% This is Q_n in lower bound on Pe,opt


function p = nctcdf0(x,nu,delta)
%NCTCDF Noncentral T cumulative distribution function (cdf).
%
%   ------------------------------------------------------------------
%   --- Modified by TCP:                                           ---
%   ---    Fixed numerical instability in computation of R1/R2:    ---
%   ---       nctcdf(159, 330, 160)  = 0 	(nonsense!)        ---
%   ---      nctcdf0(159, 330, 160)  = .4267    (correct!)         ---
%   ------------------------------------------------------------------
%
%
%   P = NCTCDF(X,NU,DELTA) Returns the noncentral T cdf with NU
%   degrees of freedom and noncentrality parameter, DELTA, at the values
%   in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also NCTINV, NCTPDF, NCTRND, NCTSTAT, TCDF, CDF.

%   References:
%      [1]  Johnson, Norman, and Kotz, Samuel, "Distributions in
%      Statistics: Continuous Univariate Distributions-2", Wiley
%      1970 p. 205.
%      [2]  Evans, Merran, Hastings, Nicholas and Peacock, Brian,
%      "Statistical Distributions, Second Edition", Wiley
%      1993 pp. 147-148.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 2.14.2.10.6.1 $  $Date: 2006/07/13 16:54:07 $

if nargin <  3,
    error('stats:nctcdf:TooFewInputs','Requires three input arguments.');
end

[errorcode x nu delta] = distchck(3,x,nu,delta);

if errorcode > 0
    error('stats:nctcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
if isa(x,'single') || isa(nu,'single') || isa(delta,'single')
   p = zeros(size(x),'single');
   seps = eps('single');
else
   p = zeros(size(x));
   seps = eps;
end

% Special cases for delta==0 and x<0.
f1 = (nu <= 0);
f0 = (delta == 0 & ~f1);
fn = (x < 0 & ~f0 & ~f1);
flag1 = any(f1(:));
flag0 = any(f0(:));
flagn = any(fn(:));
if (flag1 || flag0 || flagn)
   fp = ~(f1 | f0 | fn);
   if flag1,        p(f1) = NaN; end
   if flag0,        p(f0) = tcdf(x(f0),nu(f0)); end
   if any(fp(:)),   p(fp) = nctcdf(x(fp), nu(fp), delta(fp)); end
   if flagn,        p(fn) = 1 - nctcdf(-x(fn), nu(fn), -delta(fn)); end
   return
end

kzero = find(x == 0);
kpos = 1:numel(x);
kpos(kzero(:)) = [];

%Value passed to Incomplete Beta function.
tmp = (x.^2)./(nu+x.^2);

% Set up for infinite sum.
done = 0;
dsq = delta.^2;

% Compute probability P[t<0] + P[0<t<x], starting with 1st term
p(:) = normcdf(-delta,0,1);

% Now sum a series to compute the second term
k0 = find(x~=0);
if any(k0(:))
   tmp = tmp(k0);
   nu = nu(k0);
   dsq = dsq(k0);
   signd = sign(delta(k0));
   subtotal = zeros(size(k0));

   % Start looping over term jj and higher, this should be near the
   % peak of the E part of the term (see below)
   jj = 2 * floor(dsq/2);

   % Compute an infinite sum using Johnson & Kotz eq 9, or new
   % edition eq 31.16, each term having this form:
   %      B  = betainc(tmp,(j+1)/2,nu/2);
   %      E  = (exp(0.5*j*log(0.5*delta^2) - gammaln(j/2+1)));
   %      term = E .* B;
   %
   % We'll compute betainc at the beginning, and then update using
   % recurrence formulas (Abramowitz & Stegun 26.5.16).  We'll sum the
   % series two terms at a time to make the recurrence work out.

   E1 =          exp(0.5* jj   .*log(0.5*dsq) - dsq/2 - gammaln( jj   /2+1));
   E2 = signd .* exp(0.5*(jj+1).*log(0.5*dsq) - dsq/2 - gammaln((jj+1)/2+1));
   B1 = betainc(tmp,(jj+1)/2,nu/2);
   B2 = betainc(tmp,(jj+2)/2,nu/2);
   R1 = exp(gammaln((jj+1)/2+nu/2)-gammaln((jj+3)/2)-gammaln(nu/2) + ...
            ((jj+1)/2) .* log(tmp) + (nu/2).*log(1-tmp));
   R2 = exp(gammaln((jj+2)/2+nu/2)-gammaln((jj+4)/2)-gammaln(nu/2) + ...
            ((jj+2)/2) .* log(tmp) + (nu/2).*log(1-tmp));
   E10 = E1; E20 = E2; B10 = B1; B20 = B2; R10 = R1; R20 = R2; j0 = jj;
   todo = true(size(dsq));
   while(true)
      %Probability that t lies between 0 and x (x>0)
      twoterms = E1(todo).*B1(todo) + E2(todo).*B2(todo);
      subtotal(todo) = subtotal(todo) + twoterms;
      
      if(twoterms == -Inf)
          error('Error: twoterms = -Inf!')
      end

      
      % Convergence test.
      todo(todo) = (abs(twoterms) > (abs(subtotal(todo))+seps)*seps);
      if (~any(todo))
         break;
      end
      

      % Update for next iteration
      jj = jj+2;

      E1(todo) = E1(todo) .* dsq(todo) ./ (jj(todo));
      E2(todo) = E2(todo) .* dsq(todo) ./ (jj(todo)+1);

      B1(todo) = B1(todo) - R1(todo);
      B2(todo) = B2(todo) - R2(todo);

      R1(todo) = R1(todo) .* tmp(todo) .* (jj(todo)+nu(todo)-1) ./ (jj(todo)+1);
      R2(todo) = R2(todo) .* tmp(todo) .* (jj(todo)+nu(todo)  ) ./ (jj(todo)+2);
   end

   % Go back to the peak and start looping downward as far as necessary.
   E1 = E10; E2 = E20; B1 = B10; B2 = B20; R1 = R10; R2 = R20;
   jj = j0;
   todo = (jj>0);
   while any(todo)
      E1(todo) = E1(todo) .* (jj(todo)  ) ./ dsq(todo);
      E2(todo) = E2(todo) .* (jj(todo)+1) ./ dsq(todo);

      R1(todo) = R1(todo) .* (jj(todo)+1) ./ ((jj(todo)+nu(todo)-1) .* tmp(todo));
      R2(todo) = R2(todo) .* (jj(todo)+2) ./ ((jj(todo)+nu(todo))   .* tmp(todo));

      B1(todo) = B1(todo) + R1(todo);
      B2(todo) = B2(todo) + R2(todo);

      twoterms = E1(todo).*B1(todo) + E2(todo).*B2(todo);
      subtotal(todo) = subtotal(todo) + twoterms;

      jj = jj - 2;
      todo(todo) = (abs(twoterms) > (abs(subtotal(todo))+seps)*seps) & ...
                   (jj(todo) > 0);
   end
   p(k0) = min(1, max(0, p(k0) + subtotal/2));
end

function f = logfntheta(n, theta)
%
%  function that works as follows (for going from theta to log M)
%
%   1/fn(theta)	<= M(theta) <= 1/fn(theta) / (1-tan(theta)^2/n)
%

if(theta >= pi/2)
	f = +Inf;
end

if(theta <= 0)
	f = -Inf;
end

f = gammaln(n/2 + 1) - gammaln(n/2 + 1/2) - log(n) - log(pi)/2 + ...
	(n-2)*log(sin(theta)) + log(tan(theta));
f = log2(exp(1)) * f;

