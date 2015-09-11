function plot_all(P,n_t,n_r,T,epsilon)

% Parameters
%n_vec = [20:20:800];      % Block lengths to plot
n_vec = [20 40 80 120 160 200 250 300 350 400 500 600];
%P = 10;                   % SNR (in DB)
%n_t = 2;                  % # transmit antennas
%n_r = 2;                  % # receiver antennas
%T = 2;                    % Coherence time
%epsilon = .001;           % (Average) Probability of error
real_or_complex = 'real';  % Specifies real or complex MIMO channel, set to 'real' for real, 'complex' for complex.

% Plotting
vals = [];
apx_vals = [];
for n = n_vec
    log_M = MIMO_achievability_CSIR(n,epsilon,n_t,n_r,T,P,real_or_complex);   
    vals = [vals 1/(n*T)*log_M]; % Note this is logm for n*T samples, not increasing in 1
    n
    1/(n*T)*log_M
end

plot(T*n_vec,vals)
grid on
% Adds normal approximation curve to the plot
hold on
[C,V] = capacity_and_dispersion(n_t,n_r,T,P,real_or_complex);
plot(T*n_vec,C - sqrt(V./(T*n_vec))*qfuncinv(epsilon),'r')
hold on
plot([0 T*n_vec], C*ones(length(n_vec)+1)) % +1 is just so the capacity line extends all the way to the y-axis
