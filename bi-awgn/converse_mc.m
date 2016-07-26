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
%   'full'    - (this option will be available soon) the full bound derived
%               numerically from [1].
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
%     http://arxiv.org/abs/1511.04629
%
% Please give credits to [1] if you use this code for your research.
%
% PS1: The function  is implemented in  MATLAB R2011b and extensively uses
%      the 'fmincon' solver.  More recent MATLAB versions may require some
%      modifications.
%
% PS2: The code is tested  on a fairly large range of values.  However, if
%      the code fails on a particular choice of parameters, please drop an
%      email to erseghe@dei.unipd.it and I will try to fix the bug.
%
% (c) T. Erseghe 2016
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

% override 'full' option which is not available at the moment
if strcmp(x3,'full')
    disp('Sorry, the "full" option is not available... switching to "On2" approximation')
    x3 = 'On2';
end


% optimization settings
options1 = optimset(...
    'Algorithm','interior-point',...
    'GradObj','on',...
    'GradCon','on',...
    'Display', 'off',...
    'TolFun',1e-20,...
    'TolX', 1e-20,...
    'MaxIter',1000);
options2 = optimset(...
    'Algorithm','interior-point',...
    'GradObj','on',...
    'GradCon','off',...
    'Display', 'off',...
    'TolFun',1e-20,...
    'TolX', 1e-20,...
    'MaxIter',1000);

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
        
        try
            if (n==inf)||strcmp(x3,'normal')
                % case 1: normal approximation
                y(k) = rate_normal(n,om,qe);
                
            elseif strcmp(x3,'On2')
                % case 2: O(n^-2) approximation
                % define the saddle point starting point
                [~, ~, V] = fun_HP(om,0);
                sb0 = -qe/sqrt(V*n);
                % find saddle point
                sb = fmincon(@(x)zero_obj(x),sb0,[],[],[],[],[],[],...
                             @(x)constr_om_on2_rate(x,om,pe,n),options2);
                % calculate rate value
                y(k) = -logFA_on2(sb,om,n)/log(2);
                
            elseif strcmp(x3,'On3')
                % case 3: O(n^-2) approximation
                % define the saddle point starting point
                [~, ~, V] = fun_HP(om,0);
                sb0 = -qe/sqrt(V*n);
                % find saddle point
                sb = fmincon(@(x)zero_obj(x),sb0,[],[],[],[],[],[],...
                             @(x)constr_om_on3_rate(x,om,pe,n),options2);
                % calculate rate value
                y(k) = -logFA_on3(sb,om,n)/log(2);
                
            elseif strcmp(x3,'full')
                % case 4: full bound (numerical)
                
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
        elseif ~isnan(y(k-1)) % use previous value
            om0 = y(k-1);
        end
        
        try 
            if (n==inf)||strcmp(x3,'normal')
                % case 1: normal approximation
                % calculate SNR value
                y(k) = fmincon(@(x)zero_obj(x),om0,[],[],[],[],[],[],...
                               @(x)constr_om_normal(x,R,qe,n),options2);

            elseif strcmp(x3,'On2')
                % case 2: O(n^-2) approximation
                % define the saddle point starting point
                [~, ~, V] = fun_HP(om0,0);
                sb0 = -qe/sqrt(V*n);
                % check that if starting point is feasible
                [~, con] = constr_om_on2_error([sb0; om0],R,pe,n);
                % find SNR/saddle point
                x = fmincon(@(x)zero_obj(x),[sb0; om0],[],[],[],[],[],[],...
                            @(x)constr_om_on2_error(x,R,pe,n),options2);
                y(k) = x(2);

            elseif strcmp(x3,'On3')
                % case 3: O(n^-3) approximation
                % define the saddle point starting point
                [~, ~, V] = fun_HP(om0,0);
                sb0 = -qe/sqrt(V*n);
                % check if the starting point is feasible
                [~, con] = constr_om_on3_error([sb0; om0],R,pe,n);
                % find SNR/saddle point
                x = fmincon(@(x)zero_obj(x),[sb0; om0],[],[],[],[],[],[],...
                            @(x)constr_om_on3_error(x,R,pe,n),options2);
                y(k) = x(2);

            elseif strcmp(x3,'full')
                % case 4: full bound (numerical)

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

% constraints for calculating the error rate bound under the normal
% approximation
function [c ceq gradc gradceq] = constr_om_normal(x,R,qe,n)
c = [];
gradc = [];
ceq = [rate_normal(n,x,qe)-R];
gradceq = []; % too long to derive !!!

% log(FA) function under the O(n^-2) approximation
function y = logFA_on2(sb,om,n)
sa = sb+1;
[Hp, la, a2] = fun_HP(om,sb);
y = la*sa + log(Hp/2)-log(2*pi*n*sa^2*a2)/2/n;

% log(MD) function under the O(n^-2) approximation
function y = logMD_on2(sb,om,n)
[Hp, la, a2] = fun_HP(om,sb);
y = la.*sb + log(Hp) -log(2*pi*n*sb.^2.*a2)/2/n;

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
