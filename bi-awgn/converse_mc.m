function y = converse_mc(N, Pe, x2, x3, x4)

% converse_mc    : Polyanskiy-Poor-Verdu (PPV) metaconverse bound for
%                  binary transmission corrupted by AWGN (bi-AWGN channel)
%
% R = converse_mc(N, Pe, SNRdB, 'approx') returns the upper bound on rate
% for a codeword of length N, a packet error rate of Pe, and a signal to
% noise ratio of SNRdB (in dB). Multiple bound values can be obtained by
% ensuring that N, Pe, and SNRdB have the same length. The specific appro-
% ximation used for calculating the bound is specified by 'approx':
%
%   'normal'  - normal approximation, the O(n^-1) approximation available
%               from Fig.6 (3rd block) in [1] ; this is a very fast 
%               algorithm, but may provide a weak approximation.
%
%   'On2'     - (default) the O(n^-2) approximation derived from Fig.6 
%               (2nd block) in [1]; this is approximately 50 times slower  
%               than the 'normal' option, but much more reliable at low Pe.
%
%   'On3'     - the O(n^-3) approximation derived from Fig.6 (2nd block)
%               in [1]; this is twice slower than 'On2', and much more
%               unrealiable at high Pe. It can  be used to identify the
%               region where 'On2' is a reliable solution, namely the 
%               region where 'On2' and 'On3' closely match.
%
%   'full'    - the full bound derived numerically as from [2].
%               Be aware this is very slow to compute! It is also
%               numerically less reliable where On2 and On3 closely agree.
%
% SNRdB = converse_mc(N, Pe, R, 'approx','error') for a codeword of length
% N, and code rate R, the function returns the SNR (in dB) at which the
% lower bound on error rate is equal to Pe. Multiple bound values can be
% obtained by ensuring that N, Pe, and R have the same length. The specific
% approximation used for calculating the bound is specified by 'approx'
% (see above).
%
% All bounds were derived from:
% [1] T. Erseghe, "Coding in the Finite-Blocklength Regime: Bounds
%     based on Laplace Integrals and their Asymptotic Approximations",
%     IEEE Transactions on Information Theory, 62(12), pp 6854-6883, 2016
%     http://arxiv.org/abs/1511.04629
% [2] T. Erseghe, N. Laurenti, M. Zecchin, "Coding Bounds in the 
%     Finite-Block-Length Regime: an Application to Spread-Spectrum 
%     Systems Design," 
%     2019 AEIT International Annual Conference.
%
% Please give credits to [1] & [2] if you use this code for your research.
%
% PS1: The function  is implemented in  MATLAB R2017a and extensively uses
%      the 'fmincon' solver.  Different  MATLAB  versions may require some
%      modifications.
%
% PS2: The code is tested  on a fairly large range of values.  However, if
%      the code fails on a particular choice of parameters, please drop an
%      email to erseghe@dei.unipd.it and I will try to fix the bug.
%
% (c) T. Erseghe 2020
%

warning ('off','all');

if ~exist('x3') % default option
    x3 = 'On2';
end

if ~exist('x4') % default option
    x4 = 'rate';
end

N = N(:).'; % put into vector form
Pe = Pe(:).'; % put into vector form
x2 = x2(:).'; % put into vector form

% ensure that all input vectors have the same length
if (length(N)~=length(Pe))||(length(N)~=length(x2))
    error(['Input vectors must have the same length'])
end

% check if correct input values are given
switch x3
    case {'normal','On2','On3','full'}
    otherwise
        error(['Wrong input ',x3])
end
switch x4
    case {'rate','error'}
    otherwise
        error(['Wrong input ',x4])
end


% optimization settings
options1 = optimoptions('fmincon',...
    'Algorithm','sqp',...
    'SpecifyObjectiveGradient',true,...
    'SpecifyConstraintGradient',true,...
    'Display', 'off',...
    'TolFun',1e-20,...
    'TolX', 1e-20,...
    'MaxIterations',1000);
options2 = optimoptions('fmincon',...
    'Algorithm','sqp',...
    'SpecifyObjectiveGradient',false,...
    'SpecifyConstraintGradient',false,...
    'Display', 'off',...
    'TolFun',1e-20,...
    'TolX', 1e-20,...
    'MaxIterations',1000);

% pre-calculate Q^(-1)(Pe)
Qe = nan(size(Pe));
for pe = unique(Pe)
    qe = fmincon(@(x)zero_obj(x),0,[],[],[],[],[],[],...
                 @(x)constr_Qe(x,pe),options1);
    Qe(Pe==pe) = qe;
end

if strcmp(x4,'rate')
    
    % bound on rate
    for k = 1:length(N) % cycle on requests
        
        n  = N(k); % blocklength
        pe = Pe(k); % packet error probability
        qe = Qe(k); % Q^(-1)(pe)
        om = 10^(x2(k)/10); % snr value

        % define the saddle point starting point
        if ~((n==inf)||strcmp(x3,'normal'))
            if (k>2)&&(~isnan(y(k-2)))&&(~isnan(sb))
                sb0 = sb;
            else
                [~, ~, V] = fun_HP(om,0);
                sb0 = -qe/sqrt(V*n);
            end
        end
        
        try
            if (n==inf)||strcmp(x3,'normal')
                % case 1: normal approximation
                y(k) = rate_normal(n,om,qe);
                
            elseif strcmp(x3,'On2')
                % case 2: O(n^-2) approximation
                % find saddle point
                sb = fmincon(@(x)zero_obj(x),sb0,[],[],[],[],[],[],...
                             @(x)constr_om_on2_rate(x,om,pe,n),options2);
                % calculate rate value
                y(k) = -logFA_on2(sb,om,n)/log(2);
                
            elseif strcmp(x3,'On3')
                % case 3: O(n^-2) approximation
                % find saddle point
                sb = fmincon(@(x)zero_obj(x),sb0,[],[],[],[],[],[],...
                             @(x)constr_om_on3_rate(x,om,pe,n),options2);
                % calculate rate value
                y(k) = -logFA_on3(sb,om,n)/log(2);
                
            elseif strcmp(x3,'full')
                % case 4: full bound (numerical)                
                % find saddle point
                sb = fmincon(@(x)zero_obj(x),sb0,[],[],[],[],[],[],...
                             @(x)constr_om_full_rate(x,om,pe,n),options2);
                % calculate rate value
                y(k) = -logFA_full(sb,om,n)/log(2)/n;
                
            end
                    
            % clear unreliable values/correct values beyond the limits
            if y(k)<1/n
                y(k) = 0;
            elseif y(k)>1
                y(k) = 1;
            elseif ~isreal(y(k))
                y(k) = nan;
            end
            
        % something did not work !
        catch me
            y(k) = nan;
        end
        
        % display result
        disp(sprintf(...
            ['biawgnPPVbound(n=%d, mode=%s, approx=%s) '...
             'R = %.3g log10(Pe) = %.3g snr=%gdB '], n, x4, x3,...
              y(k), log10(pe), x2(k)));
        
    end % end cycle on k
    
elseif strcmp(x4,'error')
    
    % bound on error probability
    for k = 1:length(N) % cycle on requests
        
        n  = N(k); % blocklength
        pe = Pe(k); % packet error probability
        qe = Qe(k); % Q^(-1)(pe)
        R  = x2(k); % rate value
        
        % define the SNR starting point
        if k==1
            om0 = 1;
            [~, ~, V] = fun_HP(om0,0);
            sb0 = -qe/sqrt(V*n);
        elseif ~isnan(y(k-1)) % use previous value
            om0 = y(k-1);
        end
        if ~((n==inf)||strcmp(x3,'normal'))
            if (k>2)&&(~isnan(y(k-2)))&&(~isnan(x(1)))
                sb0 = x(1);
            else
                [~, ~, V] = fun_HP(om0,0);
                sb0 = -qe/sqrt(V*n);
            end
        end
            
        try 
            if (n==inf)||strcmp(x3,'normal')
                % case 1: normal approximation
                % calculate SNR value
                y(k) = fmincon(@(x)zero_obj(x),om0,[],[],[],[],[],[],...
                               @(x)constr_om_normal(x,R,qe,n),options2);

            elseif strcmp(x3,'On2')
                % case 2: O(n^-2) approximation
                % check that if starting point is feasible
                [~, con] = constr_om_on2_error([sb0; om0],R,pe,n);
                % find SNR/saddle point
                x = fmincon(@(x)zero_obj(x),[sb0; om0],[],[],[],[],[],[],...
                            @(x)constr_om_on2_error(x,R,pe,n),options2);
                y(k) = x(2);

            elseif strcmp(x3,'On3')
                % case 3: O(n^-3) approximation
                % check if the starting point is feasible
                [~, con] = constr_om_on3_error([sb0; om0],R,pe,n);
                % find SNR/saddle point
                x = fmincon(@(x)zero_obj(x),[sb0; om0],[],[],[],[],[],[],...
                            @(x)constr_om_on3_error(x,R,pe,n),options2);
                y(k) = x(2);

            elseif strcmp(x3,'full')
                % case 4: full bound (numerical)
                % check if the starting point is feasible
                [~, con] = constr_om_full_error([sb0; om0],R,pe,n);
                % find SNR/saddle point
                x = fmincon(@(x)zero_obj(x),[sb0; om0],[],[],[],[],[],[],...
                            @(x)constr_om_full_error(x,R,pe,n),options2);
                y(k) = x(2);
                
            end
        
        % something did not work !
        catch me
            y(k) = nan;
        end
        
        % display result
        disp(sprintf(...
            ['biawgnPPVbound(n=%d, mode=%s, approx=%s) '...
            'snr=%gdB log10(Pe) = %.3g R = %.3g '], n, x4, x3,...
            10*log10(y(k)), log10(pe), R));
        
    end % end cycle on k
    
    % turn into dBs
    y = 10*log10(y);
end

% ------------------- FUNCTIONS FOR OPTIMIZATION -------------------------

% zero object function template (for optimization under constraints only)
function [f gradf] = zero_obj(x)
f = [0];
gradf = [0];

% constraint for calculationg Q^(-1)(Pe)
function [c ceq gradc gradceq] = constr_Qe(x,Pe)
c = [];
gradc = [];
ceq = [0.5*erfc((x)/sqrt(2))-Pe];
gradceq = [-exp(-x.^2/2)/sqrt(2*pi)];

% calculates HP values
function [Hp, la, a2, a3, a4] = fun_HP(om,sb)
IL = om*(1-2*max(sb,0))-10*sqrt(om); % lower limit for integration
IU = om*(1+2*max(sb,0))+10*sqrt(om); % upper limit for integration
tol = 1e-12; % tolerance for integrations
h = @(x) log(1+exp(-2*x));
fHp = @(x) exp(-(x-om).^2/2/om-sb*h(x))/sqrt(2*pi*om);
if sb==0 % we simplify derivation
    Hp = 1;
else
    Hp = quadl(@(x)fHp(x),IL,IU,tol);
end
Hp1 = -quadl(@(x)fHp(x).*h(x),IL,IU,tol);
la = -Hp1./Hp;
if nargout<=2
    return
end
Hp2 = quadl(@(x)fHp(x).*h(x).^2,IL,IU,tol);
a2  = Hp2./Hp-(Hp1./Hp).^2;
if nargout<=3
    return
end
Hp3 = -quadl(@(x)fHp(x).*h(x).^3,IL,IU,tol);
a3  = (Hp3./Hp-3*Hp2.*Hp1./(Hp).^2+2*(Hp1./Hp).^3)/3;
if nargout<=4
    return
end
Hp4 = quadl(@(x)fHp(x).*h(x).^4,IL,IU,tol);
a4  = (Hp4./Hp-4*Hp3.*Hp1./(Hp).^2-3*(Hp2./Hp).^2 ...
    +12*Hp2.*(Hp1)^2./(Hp).^3-6*(Hp1./Hp).^4)/12;

% calculates rate under the normal approximation
function R = rate_normal(n,om,qe)
if (n==inf)
    [~, la] = fun_HP(om,0);
    C = 1 - la/log(2); % capacity
    R = C; % rate
else
    [~, la, V] = fun_HP(om,0);
    C = 1 - la/log(2); % capacity
    R = C - qe*log2(exp(1))*sqrt(V/n) + log2(n)/2/n; % rate
end

% log(FA) function under the O(n^-2) approximation
function y = logFA_on2(sb,om,n)
sa = sb+1;
[Hp, la, a2] = fun_HP(om,sb);
y = la*sa + log(Hp/2)-log(2*pi*n*sa^2*a2)/2/n;

% log(MD) function under the O(n^-2) approximation
function y = logMD_on2(sb,om,n)
[Hp, la, a2] = fun_HP(om,sb);
y = la.*sb + log(Hp) -log(2*pi*n*sb.^2.*a2)/2/n;

% log(FA) function under the O(n^-3) approximation
function y = logFA_on3(sb,om,n)
sa = sb+1;
[Hp, la, a2, a3, a4] = fun_HP(om,sb);
y = la*sa + log(Hp/2)-log(2*pi*n*sa^2*a2)/2/n;
g = (12*a2.*a4-15*a3.^2)./a2.^3/8 - 3*a3./(2*sa.*a2.^2)-1./(a2.*sa.^2);
y = y + log(1+g/n)/n;

% log(MD) function under the O(n^-3) approximation
function y = logMD_on3(sb,om,n)
[Hp, la, a2, a3, a4] = fun_HP(om,sb);
y = la.*sb + log(Hp) -log(2*pi*n*sb.^2.*a2)/2/n;
g = (12*a2.*a4-15*a3.^2)./a2.^3/8 - 3*a3./(2*sb.*a2.^2)-1./(a2.*sb.^2);
y = y + log(1+g/n)/n;

% log(FA) function under the full approximation
function y = logFA_full(sb,om,n)
[~, la, a2, a3] = fun_HP(om,sb);
% integration path
t0 = 1; % this is set heuristically
pb = @(t) sb + 1i*t/sqrt(a2) + a3*((t<=t0).*t.^2+(t>t0)*t0^2)/(2*a2^2);
pb1 = @(t) 1i/sqrt(a2) + a3*(t<=t0).*t/(a2^2);
% integrals
h = @(x) log(1+exp(-2*x));
I1 = 40;
H = @(s) (1/sqrt(2*pi*om))*quadl(@(x)exp((-1/(2*om))*(x-om).^2-s*h(x)),-I1,I1);
Beta = @(s) 2*(la*s+log(H(s)));
Alpha = @(s) Beta(s-1)-2*log(2)+2*la;
I = 40; % this is set heuristically
tol = 1e-12; % tolerance for integrations
y = integral(@(u)imag(exp(Alpha(pb(u)+1)*n/2).*pb1(u)./(1+pb(u)))/pi,0,I,'AbsTol',tol,'ArrayValued',true);
y = y + 1.0*(sb<-1.0) + 0.5*(sb==-1.0);
y = log(y);

% log(MD) function under the full approximation
function y = logMD_full(sb,om,n)
[~, la, a2, a3] = fun_HP(om,sb);
% integration path
t0 = 1; % this is set heuristically
pb = @(t) sb + 1i*t/sqrt(a2) + a3*((t<=t0).*t.^2+(t>t0)*t0^2)/(2*a2^2);
pb1 = @(t) 1i/sqrt(a2) + a3*(t<=t0).*t/(a2^2);
% integrals
h = @(x) log(1+exp(-2*x));
I1 = 40;
H = @(s) (1/sqrt(2*pi*om))*quadl(@(x)exp((-1/(2*om))*(x-om).^2-s*h(x)),-I1,I1);
Beta = @(s) 2*(la*s+log(H(s)));
I = 40; % this is set heuristically
tol = 1e-12; % tolerance for integrations
y = integral(@(u)imag(exp(Beta(pb(u))*n/2).*pb1(u)./(-pb(u)))/pi,0,I,'AbsTol',tol,'ArrayValued',true);
y = y + 1.0*(sb>0.0) + 0.5*(sb==0.0);
y = log(y);

% constraints for calculating the error rate bound under the normal
% approximation
function [c ceq gradc gradceq] = constr_om_normal(x,R,qe,n)
c = [];
gradc = [];
ceq = [rate_normal(n,x,qe)-R];
gradceq = []; % too long to derive !!!

% constraints for calculating the code rate bound under the O(n^-2)
% approximation
function [c ceq gradc gradceq] = constr_om_on2_rate(sb,om,Pe,n)
c = [];
gradc = [];
ceq = [logMD_on2(sb,om,n)-log(Pe)/n];
gradceq = []; % too long to derive !!!

% constraints for calculating the error rate bound under the O(n^-2)
% approximation
function [c ceq gradc gradceq] = constr_om_on2_error(x,R,Pe,n)
c = [];
gradc = [];
ceq = [logMD_on2(x(1),x(2),n)-log(Pe)/n; ...
    logFA_on2(x(1),x(2),n)+R*log(2)];
gradceq = []; % too long to derive !!!

% constraints for calculating the code rate bound under the O(n^-3)
% approximation
function [c ceq gradc gradceq] = constr_om_on3_rate(sb,om,Pe,n)
c = [];
gradc = [];
ceq = [logMD_on3(sb,om,n)-log(Pe)/n];
gradceq = []; % too long to derive !!!

% constraints for calculating the error rate bound under the O(n^-3)
% approximation
function [c ceq gradc gradceq] = constr_om_on3_error(x,R,Pe,n)
c = [];
gradc = [];
ceq = [logMD_on3(x(1),x(2),n)-log(Pe)/n; ...
    logFA_on3(x(1),x(2),n)+R*log(2)];
gradceq = []; % too long to derive !!!

% constraints for calculating the code rate bound under the full
% approximation
function [c ceq gradc gradceq] = constr_om_full_rate(sb,om,Pe,n)
c = [];
gradc = [];
ceq = [logMD_full(sb,om,n)-log(Pe)];
gradceq = []; % too long to derive !!!

% constraints for calculating the error rate bound under the full
% approximation
function [c ceq gradc gradceq] = constr_om_full_error(x,R,Pe,n)
c = [];
gradc = [];
ceq = [logMD_full(x(1),x(2),n)-log(Pe); ...
    logFA_full(x(1),x(2),n)+n*R*log(2)];
gradceq = []; % too long to derive !!!
