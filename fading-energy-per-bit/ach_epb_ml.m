%Info_bits = ach_epb_ml(EE, epsil)
%
%This function computes the achievability bound on the maixmal number of
%information bits that can be transmitted over a Rayleigh-fading channel
%(no CSI) with a fixed energy budget "EE" and error probability "epsil".
%Here, "EE" is a vector, and "epsil" is a scalar. Note that the entries in
%"EE" must be arranged in ascending order.     

function info_bits=ach_epb_ml(EE,epsil)


%normal approximation for k^*(E,error) 
info_bits_na=EE/log(2) - (12^(-1/3)+(2/3)^(1/3))*(qfuncinv(epsil)/log(2))^(2/3)*(log2(EE)).^(1/3).*(EE).^(2/3);

info_bits=[];

num_samples=1000; %number of k_values

%determine the optimal number of nonzero entries using normal approximation
NN_na = (1.5*qfuncinv(epsil)*EE/log(2)./log2(EE)).^(2/3);

exceeding_indicator=0; %1 means that E is too big to do exact computations.

loop=10000;

for index_E = 1:length(EE)
    
    E = EE(index_E);
    disp(['E=',num2str(E)]);
    tic
    
    N_na = NN_na(index_E); %optimal N from normal approximation
    
    N_start = floor(N_na/5);
    N_end = floor(N_na*5);
    
    
    step = (log(N_end) -log(N_start))/500;
    
    NN = [round(exp(log(N_start):step: log(N_end))), N_end];
    NN=unique(NN);%remove multipicity
    
    if exceeding_indicator==0
        
        k_start = max(floor(info_bits_na(index_E)/2),1);
        
        k_end = info_bits_na(index_E)*2;
        
        step_k = (log(k_end)-log(k_start))/num_samples;
        
        kkk_ary = exp(log(k_start):step_k:log(k_end));
        
        index_kk =find(kkk_ary<=1);
        kkk_ary(index_kk)=[];
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
    Eb=log10(E/max(opt_k_ary) )*10;
    disp(['k^*(E,error)=', num2str(max(opt_k_ary)),';    Eb=',num2str(Eb)]);
    toc
end
%semilogx(info_bits, log10(EE./info_bits)*10)