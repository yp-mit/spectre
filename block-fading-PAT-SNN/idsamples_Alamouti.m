function i_s = idsamples_Alamouti(Mr, L, np, nc, rho, N_MC)

if nargin < 6
    N_MC = 1e6;                 % cap the number of iterations  
end

Mt = 2;                         % number of tx antennas

% Other parameters
s = 1;                          % Gallager exponent (to be optimised)
sigmaEst = 1/sqrt(rho*np/Mt);   % estimation error standard deviation
nd = nc - np;                   % data symbols per coherence time

% Generate samples
i_s = zeros(L, N_MC);
for n_mc = 1:N_MC
    for l = 1:L
        % Generate data and channel for the coherence block
        U_l = sqrt(rho/Mt) * qpsksample(2,nd/2); % symbol matrix
        u_l = reshape(U_l,nd,1);
        H_l = randcn(Mt,Mr);
        REst_l = zeros(Mr,nd);
        gEst_l = zeros(Mr,1);
        
        % For each rx antenna
        for m = 1:Mr
            % create equivalent matrix V_ml
            h_ml = H_l(:,m);
            V_ml = [h_ml(1)         h_ml(2); 
                    conj(h_ml(2))   -conj(h_ml(1))];
            % true signal received
            N_ml = randcn(2,nd/2);
            Y_ml = V_ml * U_l + N_ml;
            % estimated channel and equivalent matrix Vhat_ml
            hEst_ml = h_ml + sigmaEst * randcn(2,1);
            VEst_ml = [hEst_ml(1)           hEst_ml(2); 
                       conj(hEst_ml(2))     -conj(hEst_ml(1))];
            % rx signal for antenna m
            REst_ml = VEst_ml'/norm(hEst_ml) * Y_ml;
            rEst_ml = reshape(REst_ml,1,[]);
            % 
            REst_l(m,:) = rEst_ml;
            gEst_l(m) = norm(hEst_ml);
        end
        
        ytilde_l = gEst_l.'/norm(gEst_l) * REst_l;
        ytilde_l = ytilde_l(:);

        % Compute the information density for the block
        aux = 0;
        for t = 1:nd
            aux = aux + idgallagerqpsknn( u_l(t), ytilde_l(t), ...
                norm(gEst_l), s, rho/Mt);
        end
        i_s( l, n_mc ) = aux;
    end
end

i_s = sum(i_s,1); % row vector    
i_s = i_s(:);   % column vector
