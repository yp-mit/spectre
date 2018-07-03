function epsilon = QPSK_Alamouti_eps(Mr, k, L, nc, fp, SNR_dB, spa)
% QPSK_ALAMOUTI_EPS
%
% Inputs
%   Mr:     no of rx antennas
%   k:      number of information bits
%   L:      no of coherence blocks
%   nc:     length of coherence block
%   fp:     fraction of block for training
%   SNR_dB: SNR in dB
%   spa:    flag for saddlepoint approximation 
%   

%% Parameters and initialisation

% Simulation parameters
RELSPREAD_MAX = 0.40;           % max relative spread of outcomes

% Parallelisation
ncores = feature('numcores');   % number of cores on current machine
nworkers = ncores;              % use all cores for the algorithm
poolobj = gcp('nocreate');      % if no pool do not create new one
if isempty(poolobj)
    parpool(nworkers); % initialise the pool 
end

% System parameters
n = L * nc;             % blocklength
R = k/n;                % rate (bits)
np = floor(fp*nc);      % number of pilots

nd = nc - np;
if mod(nd,2) == 1
    np = np + 1;
end

%% Compute error probability
    
rho_dB = SNR_dB;
rho = 10^(rho_dB/10);

eps_trial = eps_RCUs_Alamouti_stable(R, Mr, L, np, nc, rho, spa, RELSPREAD_MAX);
epsilon = median( eps_trial );
