%[info_bits,opt_N] = normal_approx_nocsi(EE, epsil)
%
%This function computes the normal approximation for the maixmal number of
%information bits that can be transmitted over a Rayleigh-fading channel
%(no CSI) with a fixed energy budget "EE" and error probability "epsil".
%Here, "EE" is a vector, and "epsil" is a scalar. The function also returns
%the optimal number of channel uses "opt_N".

function [info_bits, opt_N]=normal_approx_nocsi(EE,epsil)

info_bits=[];
opt_N=[];
for E = EE
    
    loop=10^5;
    
    if E<loop
        NN=1:loop;
    else
        NN=round(exp(0:(log(E)/loop):log(E) ));
        NN=unique(NN);%remove multipicity
    end
    
    logM = E - log(1+E./NN).*NN - sqrt(E^2./NN)*qfuncinv(epsil);
    
    info_bits=[info_bits, max(logM)/log(2)];
    
    optn = find(logM==max(logM));
    
    if length(optn) ==1 
        opt_N=[opt_N, optn];
    else
        opt_N=[opt_N,optn(1)];
    end
    
end
