function i_s = idsamples(Mr, L, np, nc, rho, N_MC )

if nargin < 6
    N_MC = 1e6;             % cap the number of iterations 
end

% Other parameters
s = 1;                      % Gallager exponent (to be optimised)
sigmaEst = 1/sqrt(rho*np);  % estimation error standard deviation
nd = nc - np;               % data symbols per coherence time

% Generate samples
i_s = zeros(L,N_MC);
for n_mc = 1:N_MC
    for l=1:L
        % Generate data and channel for the coherence block
        x_l = sqrt(rho) * qpsksample(nd); % nd iid QPSK samples
        h_l = randcn(Mr,1);
        hEst_l = h_l + sigmaEst * randcn(Mr,1);
        n_l = randcn(Mr,nd);
        y_l = h_l * x_l.' + n_l;

        htildeEst_l = norm( hEst_l );
        ytilde_l = (hEst_l/norm(hEst_l))' * y_l;
        ytilde_l = ytilde_l(:);

        % Compute the information density for the block
        aux = 0;
        for t = 1:nd
            aux = aux + idgallagerqpsknn( x_l(t), ytilde_l(t), ...
                htildeEst_l, s, rho);
        end
        i_s( l, n_mc ) = aux;
    end
end

i_s = sum(i_s,1); % row vector    
i_s = i_s(:);   % column vector
    

