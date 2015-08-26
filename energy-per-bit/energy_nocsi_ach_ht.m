%
%Info_bits = ach_epb_ht(EE, epsil)
%This function computes the kappa-beta version of the achievability bound
%on the maixmal number of information bits that can be transmitted over a
%Rayleigh-fading channel with a fixed energy budget "EE" and error
%probability "epsil". Here, EE is a vector, and epsil is a scalar.   

function info_bits = energy_nocsi_ach_ht(EE, epsil)

info_bits=[];
loop=100000;

for E=EE
    %tic  
    N_max = E;  % the maximal number of nonzero entries in a codeword
    
    if N_max>1000 
        
        step = log(N_max)/1000;
        
        NN = [round(exp(0:step: log(N_max)-step)), N_max];
        
    else
        NN=1:1:N_max;
    end
    
    info_bits_tem=[];
    
    for tau = [epsil/20 epsil/15 epsil/10 epsil/5 epsil/4 epsil/3 epsil/2 epsil/1.5] %optimize over tau      
        for N = NN %optimize over N, the number of nonzero entries in a codeword
            x_th = gammaincinv(epsil-tau,N);
            k_temp = log(tau)-log(gammainc( (1+E/N)*x_th, N, 'scaledupper')) + (1+E/N)*x_th + gammaln(N+1)-N *log(x_th) - N*log(1+E/N);
            info_bits_tem = [info_bits_tem, k_temp];
        end
    end
    info_bit =  max(info_bits_tem)/log(2);
    if(isempty(info_bit))
	    info_bit = 0;
    end
    if info_bit <100 
        info_bit = log2(2^info_bit + 1);% M-1 >= tau/beta.
    end
    info_bits =[info_bits, info_bit]; % M >= tau/beta.
    %toc
end
