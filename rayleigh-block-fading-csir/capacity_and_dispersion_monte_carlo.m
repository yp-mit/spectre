function [C,V] = capacity_and_dispersion(n_t,n_r,T,P,real_or_complex)

% This function computes the capacity and dispersion of the MIMO block
% fading channel.  Rayleign fading is assumed.
%
% Parameters:
%   n_t - Number of transmit antennas
%   n_r - Number of receiver antennas
%   T   - Coherence time
%   P   - Power (note dB's = 10^(P/10))
%
% The expressions are
%   C = E[log det( I_{n_r} + P/n_t HH^T)]
%   V = Var[i(X;Y,H)|X]
%

%% Capacity
% Sampling method: note HH^T ~ 1/n_t*Wishart(I_{n_r},n_t), which matlab can
% sample from.

complex = 0;
if (nargin == 5)
    if (strcmp(real_or_complex,'complex') == 1)
        complex = 1;
    end
end

samples = 10^4;

V = eye(n_r);
for i=1:samples
    % use wishart distribution, i.e. hh' where h = 1/sqrt(n_t)*randn(n_r,n_t)
    if (complex == 0)
        hh_term = 1/n_t*wishrnd(V,n_t);
        s(i) = 1/2*log2(det(eye(n_r) + P/n_t*hh_term ));
    else
        h = 1/sqrt(n_t)*sqrt(.5)*(randn(n_r,n_t) + 1i*randn(n_r,n_t));
        %h = 1/sqrt(n_t)*(randn(n_r,n_t) + 1i*randn(n_r,n_t));
        s(i) = real(log2(det(eye(n_r) + P/n_t*h*h')));
    end
end
C = mean(s);


%% Dispersion: 
% Computing 1/T*Var(i(X;Y,H)|X) = 1/T*E_X[(E[(i(x;Y,H)] - E[i(x;Y,H))^2]]
samples1 = 500;
samples2 = 500;

if (complex == 0)
exp_x_samples = zeros(samples1,1);
for i=1:samples1 % run time = O(samples1 * samples2 * T)
    x = sqrt(P/n_t)*randn(n_t,T);
    info_dens_x_samples = zeros(samples2,1);
    
    % Compute mu_x = E[i(x;Y,H)] for this x
    for j=1:samples2
        info_dens_x = zeros(1,T);
        h = 1/sqrt(n_t)*randn(n_r,n_t);
        K = eye(n_r) + P/n_t*h*h'; 
        for k=1:T
            z = randn(n_r,1);
            y = h*x(:,k) + z;
            info_dens_x(k) = 1/2 * ( log2(det(K)) - (norm(y-h*x(:,k))^2 + y'*inv(K)*y)*log2(exp(1)) );
        end
        info_dens_x_samples(j) = sum(info_dens_x);
    end
    mu_x = mean(info_dens_x_samples);
    exp_x = mean( (info_dens_x_samples - mu_x).^2 );
    exp_x_samples(i) = exp_x;
end
V = 1/T * mean(exp_x_samples);

else
    
exp_x_samples = zeros(samples1,1);
for i=1:samples1 % run time = O(samples1 * samples2 * T)
    x = sqrt(P/n_t)*sqrt(.5)*(randn(n_t,T) + 1i*randn(n_t,T)); % iid Gaussian x
    info_dens_x_samples = zeros(samples2,1);
    
    % Compute mu_x = E[i(x;Y,H)] for this x
    for j=1:samples2
        info_dens_x = zeros(1,T);
        h = 1/sqrt(n_t)*sqrt(.5)*(randn(n_r,n_t) + 1i*randn(n_r,n_t));
        K = eye(n_r) + P/n_t*h*h'; 
        for k=1:T
            z = sqrt(.5)*(randn(n_r,1) + 1i*randn(n_r,1));
            y = h*x(:,k) + z;
            info_dens_x(k) = log2(det(K)) + (-norm(y-h*x(:,k))^2 + y'*inv(K)*y)*log2(exp(1));
        end
        info_dens_x_samples(j) = real(sum(info_dens_x)); % each entry is a sample of i(X;Y,H) for one T block
    end
    mu_x = mean(info_dens_x_samples);
    exp_x = mean( (info_dens_x_samples - mu_x).^2 );
    exp_x_samples(i) = exp_x;
end
V = 1/T * mean(exp_x_samples);

end
    
 

