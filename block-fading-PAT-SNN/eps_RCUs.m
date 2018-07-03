function eps = eps_RCUs( i_s, n, r )
eps = mean( exp( - max(0, i_s - log(2^(n*r)-1) ) ) ); % avg pr err 
% eps = mean( exp( - max(0, i_s - log(2^(n*r)-1)-log(2)) ) ); % max pr err
end