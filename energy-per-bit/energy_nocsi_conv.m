%info_bits=function conv_epb_nocsi(EE, epsil);
%
%This function computes the converse bound on the maixmal number of
%information bits that can be transmitted over a Rayleigh-fading channel
%(no CSI) with a fixed energy budget "EE" and error probability
%"epsil". Here, "EE" is a vector, and "epsil" is a scalar.


function info_bits=energy_nocsi_conv(EE, epsil)

info_bits=[];

[info_bits_na, NN_na]=energy_nocsi_normapx(EE,epsil);

if epsil<0.1
    target_error = epsil*2;
else
    target_error=epsil*1.1;
end


eps_small =10^(-8);

for ij = 1:length(EE)

    
    
    E=EE(ij);
    
    disp(['E=',num2str(E)]);
    
    N_start = NN_na(ij)/5;
    
    N_end=NN_na(ij)*5;
    
    step1 = (log(N_end) -log(N_start))/10000;
    
    NN1 = [round(exp(log(N_start):step1: log(N_end))), N_end];
    
    NN1=unique(NN1);%remove multipicity
    
    
    
    %As a first step, we want to guess a good eta. We do this by assuming
    %that the opt codeword is [x,...,x,0,...]
    
    yy=gammaincinv(target_error,NN1);
    
    eta_arr = E*yy./NN1 - NN1.*log(1+E./NN1);
    
    eta=max(eta_arr);
    
    CDF_array=[];
    
    %we first consder the scenario that the codeword contains three
    %distinct values
    %
    
    step2 = (log(N_end) -log(N_start))/500;
    
    NN = [round(exp(log(N_start):step2: log(N_end))), N_end];
    
    NN=unique(NN);
    
    tic

    for N=NN
        
        x1_min  = E/(1+N); %in this case x1= x2, this is to make sure that x1>=x2>=x3.
        
        step_x1 = E/100;
        
        for x1 = x1_min : step_x1: E-step_x1
            
            x3_max = (E-x1)/(N+1); % the max of x3 is equal to x2
            
            for x3 = 0: x3_max/10: x3_max*9/10 % we only take 10 values for x3_max
                
                x2= (E-x1-x3)/N;
                
                if x2>=x1-eps_small | x3>=x2 -eps_small | x3<= eps_small
                    continue
                end
                
                R_prob = eta + log(1+x1) + N*log(1+x2) + log(1+x3);

                CDF = gammainc(R_prob/x2,N) - x1/(x1-x3)* exp(-R_prob/x1 - N*log(1 -x2/x1) )* gammainc(R_prob*(1/x2-1/x1),N);
                arg_int=0: R_prob/10000:R_prob;
                
                pdf_GN = exp((N-1)*log(arg_int) - arg_int/x2 - N*log(x2) - gammaln(N) );
                
                int2= sum( exp(-(R_prob-arg_int)/x3).* pdf_GN  )*(R_prob/10000);
                
                CDF = CDF + x3/(x1-x3)* int2;
                
                CDF_array=[CDF_array,CDF];
         
                
            end % end x3
            
        end %end x1
        
        
    end %end N
    
    step3 = (log(N_end) -log(N_start))/50;
    
    NN2 = [round(exp(log(N_start):step3: log(N_end))), N_end];
    
    NN2=unique(NN2);
  
   
    for N1=NN2
        
        for N2=NN2
            step_x1= (E/N1 -E/(N1+N2))/10;
            for x1= E/(N1+N2): step_x1: E/N1
                x2 = (E-x1*N1)/N2;
                if x2<eps_small | x2>=x1-eps_small
                    continue;
                end
                R_prob = eta + N1*log(1+x1) + N2*log(1+x2);
                arg_int=0: R_prob/10000:R_prob;
                PDF_G2 = exp((N2-1)*log(arg_int) - arg_int/x2 - N2*log(x2) - gammaln(N2));
                CDF_G1 = gammainc((R_prob-arg_int)/x1, N1);
                
                CDF = sum(PDF_G2.*CDF_G1)*(R_prob/10000);
                CDF_array=[CDF_array,CDF];
   
                
            end
        end
        
    end
    
    toc
    
    
    indexx = find(CDF_array<0);
    CDF_array(indexx)=[];

    CDF_min = min(min(CDF_array), target_error);
    
    info_k =  (eta - log(CDF_min -epsil))/log(2);
    info_bits = [info_bits, info_k];
end
