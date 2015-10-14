function [rate_vec,C,V] = MIMO_csir_normapx(epsilon,n_t,n_r,T,P,real_or_complex,n_vec_in)

% This function computes the normal approximation to the MIMO (real or
% complex) block fading channel
%
% Returns rate_vals - a vector with entries 1/(n*T)*log M(n,epsilon,P), for
% n from n_min to n_max

% n_min = 20;
% n_max = 800;
% n_increment = 20;
% n_vec = [n_min:n_increment:n_max];
n_vec = [20 40 80 120 160 200 250 300 350 400 500 600];

complex = 1;
if (nargin >= 6) && (~isempty(real_or_complex))
    if (strcmp(real_or_complex,'real') == 1)
        complex = 0;
    end
end

if (nargin >= 7) && (~isempty(n_vec_in))
    n_vec = n_vec_in;
end

disp('Computing Integrals for Capacity and Dispersion...')
if (complex == 0)
    [C,V] = capacity_and_dispersion_mc(n_t,n_r,T,P,'real');
else
    [C,V] = capacity_and_dispersion(n_t,n_r,T,P)
end
disp('Done integrating') 

rate_vec = C - sqrt(V./n_vec/T)*qfuncinv(epsilon);
%plot(n_vec,logM_vals);
