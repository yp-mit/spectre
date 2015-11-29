%Info_bits = energy_nocsi_ach_ml(EE, epsil)
%
%This function computes the achievability bound on the maixmal number of
%information bits that can be transmitted over a Rayleigh-fading channel
%(no CSI) with a fixed energy budget "EE" and error probability "epsil".
%Here, "EE" is a vector, and "epsil" is a scalar. Note that the entries in
%"EE" must be arranged in ascending order.     

function info_bits=energy_nocsi_ach_ml(EE,epsil)

%normal approximation for k^*(E,error) 
[info_bits_na, NN_na] = energy_nocsi_normapx(EE,epsil);

info_bits=[];

num_samples=1000; %number of k_values

exceeding_indicator=0; %1 means that E is too big to do exact computations.

loop=10000;

for index_E = 1:length(EE)
    
    E = EE(index_E);
    disp(['E=',num2str(E)]);
    
    tic
    
    if info_bits_na(index_E) <=0 
%        N_na = E; % in this case, E is small, and the normal approximation may not be accurate
    disp(['This value of E is too small, energy_nocsi_ach_ml(', num2str(E),',',num2str(epsil),')=0']);
    
    info_bits=[info_bits, 0];
    
    continue;
    else
        N_na = NN_na(index_E); %optimal N from normal approximation    
    end
    
    
    N_start = floor(N_na/5);
    N_end = floor(N_na*5);
    
    
    step = (log(N_end) -log(N_start))/500;
    
    NN = [round(exp(log(N_start):step: log(N_end))), N_end];
    NN=unique(NN);%remove multipicity
    
    if exceeding_indicator==0
        
        k_start = max(floor(info_bits_na(index_E)/2),1);
        
        k_end = info_bits_na(index_E)*5;
        
        step_k = (log(k_end)-log(k_start))/num_samples;
        
        kkk_ary = exp(log(k_start):step_k:log(k_end));
        
        index_kk =find(kkk_ary<=1);
        kkk_ary(index_kk)=[];
        
        
        kkk_ary = [kkk_ary, ceil(k_end):1: 1030]; % this is to make sure that the "while" loop ends before kkk_ary was exhausted 
        
    end
    
    opt_k_ary=[];
    
    for index_N = 1:length(NN) %optimize over N
        
        N=NN(index_N);        
        ij=1;
        while exceeding_indicator==0 
            
            if kkk_ary(ij) <= 1020 %2^(1020) is probably the largest real number that a 64-bit version of MATLAB can handle
                M = 2.^(kkk_ary(ij));
                x_th = gammaincinv(1/(M-1), N, 'upper')/(1+E/N); %(M-1)*prob[\bar{G}_N> (1+E/N)x_th]=1, where $\bar{G}_N~Gamma(N,1)$
                
                step_z = 0.01;
                zz= 1:step_z:100;
                
                integrand = gammainc(x_th*zz, N) .* exp( (N-1).*log(zz) - (zz-1)*x_th*(1+E/N) );
                
                error_M = N * sum(integrand)*step_z/gammainc(x_th*(1+E/N),N,'scaledupper');

                if error_M>epsil
                    if ij==1
                        opt_k_ary = [opt_k_ary,0];
                    else
                        opt_k_ary=[opt_k_ary, kkk_ary(ij-1)];
                    end
                    break;
                end
                
            else
                % in this case it is difficult to do exact computation. Instead, we should bound gammainc.
                exceeding_indicator=1;
            end
            ij=ij+1;
        end
        
        
        if exceeding_indicator == 1
            xx= N/loop: N/loop: N;
            
            index =sum( (gammainc(xx,N) + 4*exp( - xx + (N-1)*log(xx) -gammaln(N))/(2+E/N)) < epsil );
            
            if index ==0
                opt_k_ary = [opt_k_ary,0];
            else
                x_th=xx(index);
                
                info_nats = x_th*(1+E/N)-(N-1)*log(x_th)-(N-1)*log(1+E/N)+gammaln(N)-log(2);
                opt_k_ary = [opt_k_ary, info_nats./log(2)];                
            end
        end
        
        if index_N>4
            delta_k = opt_k_ary(end-2:end) - opt_k_ary(end-3:end-1);
            if sum(delta_k <0) == 3
                break;
            end
        end
    end%end N
    
    info_bits=[info_bits, max(opt_k_ary)]; 
    toc
end

