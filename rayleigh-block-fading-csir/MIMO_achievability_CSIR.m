function log_M = MIMO_achievability_CSIR(n,epsilon,n_t,n_r,T,P,real_or_complex)

% This computes beta_{alpha}(P_{XY},P_XQ_Y) when P_X is the Haar measure
% on the input space and Q_Y is the CAOD for the MIMMO block fading channel
%
%  P_X = X^n * sqrt(nTP) / ||X^n||_F, X_{ij} ~ N(0,P/n_t)
%  Q_{Y_i} = P_H * N(0, I_{n_r} + P/n_t HH^T),   i = 1,...,nT
%
% We assume Rayleigh fading.  We use the average probabily of error
% achievability bound listed in the documentation.  A short summary of the
% code is: we compute (via monte carlo)
%
%   log(gamma) s.t. P_{XY}[log P_{XY}/P_XQ_Y >= log(gamma)] = alpha
%   output: P_XP_Y[log P_{XY}/P_XQ_Y >= log(gamma)]
%
% Note that we only need to sample log P_{XY}/P_XQ_Y under P_{XY}, since
%
%   P_XQ_Y[sum i(x;y,h) > log_gamma] 
%       = E[e^{-log(P_{XY}/P_XQ_Y)1(\sum i(x;y,h) > log_gamma)]
%
% Where the expecation is w.r.t. P_{XY}

%%%%% MONTE CARLO %%%%%

% Sample from log(P_XY / P_XQ_Y) under P_XY
iter = 10^3;  % # of monte carlo samples

complex = 0;  % Set this to 1 for complex MIMO, otherwise 0 for real case
if (nargin == 7)
    if (strcmp(real_or_complex,'complex') == 1)
        complex = 1;
    end
end


if (complex == 0)
% Sample x uniform on {x \in R^{n_t x T*n} : ||x||_F^2 = nTP}
x = zeros(n_t,T,n,iter);
for i=1:iter
    x_gauss = randn(n_t,n*T);
    for j = 1:n
        x_mat = x_gauss * sqrt(n*T*P) / norm(x_gauss,'fro');
        x(:,:,j,i) = x_mat(:,(j-1)*T+1:j*T);
    end
end
%x = sqrt(P/n_t)*randn(n_t,T,n,iter); % iid Gaussian P_x
h = 1/sqrt(n_t)*randn(n_r,n_t,n,iter);  % Fading coefficients
z = randn(n_r,T,n,iter);                % Additive noise
y = zeros(n_r,T,n,iter);                % Channel output
K = zeros(n_r,n_r,n,iter);              % Output Covariance matrix
info_density_samples = zeros(n*T,iter); % Matrix of info density samples

else
% Same as above, except complex case    
x = zeros(n_t,T,n,iter);
for i=1:iter
    x_gauss = sqrt(.5*P/n_t)*(randn(n_t,n*T) + 1i*randn(n_t,n*T));
    for j = 1:n
        x_mat = x_gauss * sqrt(n*T*P) / norm(x_gauss,'fro');
        x(:,:,j,i) = x_mat(:,(j-1)*T+1:j*T);
    end
end
%x = sqrt(P/n_t)*randn(n_t,T,n,iter); % iid Gaussian P_x
h = 1/sqrt(n_t)*sqrt(.5)*(randn(n_r,n_t,n,iter) + 1i*randn(n_r,n_t,n,iter));  % Fading coefficients
z = sqrt(.5)*(randn(n_r,T,n,iter) + 1i*randn(n_r,T,n,iter));                  % Additive noise
y = zeros(n_r,T,n,iter);                % Channel output
K = zeros(n_r,n_r,n,iter);              % Output Covariance matrix
info_density_samples = zeros(n*T,iter); % Matrix of info density samples

end % end complex vs real logic

for i=1:n
    for j=1:iter
        y(:,:,i,j) = h(:,:,i,j)*x(:,:,i,j) + z(:,:,i,j);
        K(:,:,i,j) = eye(n_r) + P/n_t*h(:,:,i,j)*(h(:,:,i,j)');
        for k=1:T
            % matrix of samples, index by nT by iter, summing over rows
            % gives sample of sum i(x^n;y^n,h^n)
            if (complex == 0)
                info_density_samples((i-1)*T+k,j) = 1/2 * ( log2(det(K(:,:,i,j))) - norm(z(:,k,i,j),'fro')^2 + y(:,k,i,j)'*inv(K(:,:,i,j))*y(:,k,i,j) );
            else % Only the factor of 2 changes in the complex case, and we take real part (the expression is real anyway, but matlab adds "+ 1i*0.0000"
                info_density_samples((i-1)*T+k,j) = real(log2(det(K(:,:,i,j))) - norm(z(:,k,i,j),'fro')^2 + y(:,k,i,j)'*inv(K(:,:,i,j))*y(:,k,i,j));
            end
        end
    end
end

n_dens_samples = sort(sum(info_density_samples),'descend'); % Sorted vector of iter samples of sum i(x;y,h)

% Find alpha-quantile of the distribution of the information density
tau = epsilon / 2;
alpha = 1 - epsilon + tau;
if (ceil(iter*alpha) == 0)    
    log_gamma = n_dens_samples(1);
else
    log_gamma = n_dens_samples(ceil(iter*alpha));
end

% Div is the divergence between P_Y induced by the uniform (Haar) input measure
% on the space {x^n \in R^{n_t x nT}: \|x^n\|_F^2 = nTP}, and Q_Y, the
% ouput measure for the i.i.d. Gaussian input measure
if (gamma((n*T*n_t+1)/2) ~= Inf)
    div = 2*P*(n*T - sqrt(2*n*T/n_t)*gamma((n*T*n_t+1)/2)/gamma(n*T*n_t/2));
else
    div = P/(2*n_t);
end
% div = P/(2*n_t); % A reasonably approximation to the above divergence
lambda = 4;  % Renyi paramter
const = lambda/(lambda-1)*log2(tau) - 1/(lambda-1)*log2(2^((lambda-1)*lambda*div) - (1-tau)^lambda); % Lower bound on Beta_{\tau}(P_Y,Q_Y)

% Compute log(beta_{alpha}(P_{XY},P_XQ_Y)) while keeping the exponential not decaying as e^{-n}
log_beta = -log_gamma + log2(mean(2.^(-(n_dens_samples - log_gamma)).*(n_dens_samples >= log_gamma)));

log_M = const - log_beta;
