function eps_trial = eps_RCUs_Alamouti_stable(R, Mr, L, np, nc, rho, spa,...
                RELSPREAD_MAX, N_MC_MIN, N_MC_MAX, N_REP_TRIALS, d_verbose)
%EPS_RCUS_ALAMOUTI_STABLE returns the error probability computed via RCUs 
%in a system using Alamouti coding (Mt=2 tx antennas) with parameters: 
%R (rate), Mr (number of rx antennas), L (number of blocks), 
%np (number of pilots), nc (size of coherence block), rho (SNR linear).

	% Default values for parameters that are not passed
    if nargin < 12
        d_verbose = 2;      % set verbosity
    end
    
    if nargin < 11
        ncores = feature('numcores');   % number of cores on current machine
        nworkers = ncores;              % use ncores-1 in parallel 
        poolobj = gcp('nocreate');      % if no pool do not create new one
        if isempty(poolobj)
            parpool(nworkers);          % initialise the pool 
        end
        % choose N_REP_TRIALS a multiple of no. of workers for high
        % efficiency
        N_REP_TRIALS = 3 * poolobj.NumWorkers; 
    end
    
    if nargin < 10
        N_MC_MAX = 2^19;    % maximum number of Monte Carlo iterations
    end
    
    if nargin < 9
        N_MC_MIN = 2^10;    % initial number of Monte Carlo iterations
    end
    
    if nargin < 8
        RELSPREAD_MAX = 0.50;   % max relative spread of outcomes
    end
    
    if nargin < 7
        spa = 0;
    end
    % Algorithm parameters
    
    N_MC = min([N_MC_MIN N_MC_MAX]);	% initial number of Monte Carlo trials
    
    % Algorithm
    n = nc * L;
    b_iter = 1;
    
    if d_verbose >= 1
        fprintf([' Estimating Pe with ' ...
            '%d trials...\n'],N_MC);
    end
    
    while b_iter == 1
        
        eps_trial = zeros(N_REP_TRIALS,1);
        
        parfor i_trial = 1:N_REP_TRIALS
            if spa == 0 
                % Generate information density samples
                i_s = idsamples_Alamouti(Mr, L, np, nc, rho, N_MC);

                % Compute error probability
                eps_trial(i_trial) = eps_RCUs(i_s, n, R);
            else %saddlepoint approximation
                % Generate information density samples with L = 1
                i_s = idsamples_Alamouti(Mr, 1, np, nc, rho, N_MC);
                
                % Compute error probability using saddlepoint approximation
                eps_trial(i_trial) = eps_RCUs_SPA(i_s, n, L, R);
            end

        end
                
        % Check spread of data
        relspread_eps = relspread(eps_trial); % relative spread
        
        if d_verbose >= 2
            fprintf('  %.3e', eps_trial);
            fprintf('\n  Relative spread: %.3e\n', relspread_eps);
        end 
                
        % Stopping conditions
        if relspread_eps < RELSPREAD_MAX % accuracy reached
            b_iter = 0;
            if d_verbose >= 1
                fprintf(' Estimation accuracy reached with %d trials.\n',N_MC);
            end
        else
            if 2 * N_MC <= N_MC_MAX
                N_MC = 2 * N_MC;
                if d_verbose >= 1
                    fprintf(['  Accuracy not reached: '...
                        'trying estimation of Pe with ' ...
                        '%d trials...\n'], N_MC);
                end
            else
                b_iter = 0;
                if d_verbose >= 1
                    warning('::Relative accuracy target not reached');
                end
            end
        end
    end
    
end % of function