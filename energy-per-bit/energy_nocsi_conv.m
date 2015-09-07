%info_bits=function energy_nocsi_conv(EE, epsil);
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


eps_small =10^(-10);

num_loop_NN=500;
num_loop_x1=100;
num_loop_x3=10;

num_loop_NN2=50;
num_loop_x12=10;

for ij = 1:length(EE)
    
    E=EE(ij);
    
    disp(['E=',num2str(E)]);
    
    if info_bits_na(ij)<0
        NN_na(ij) = EE(ij);
    end
    
    N_start = max(floor(NN_na(ij)/5),1);
    N_end = ceil(NN_na(ij)*5);
    step1 = (log(N_end) -log(N_start))/10000;
    NN1 = [round(exp(log(N_start):step1: log(N_end))), N_end];
    NN1=unique(NN1);%remove multipicity
    
    %As a first step, we want to guess a good eta. We do this by assuming
    %that the opt codeword is [x,...,x,0,...]
    yy=gammaincinv(target_error,NN1);
    eta_arr = E*yy./NN1 - NN1.*log(1+E./NN1);
    eta=max(eta_arr);
   
    %we first consider the scenario that the codeword contains three
    %distinct values 

    step2 = (log(N_end) -log(N_start))/num_loop_NN;
    NN = [round(exp(log(N_start):step2: log(N_end))), N_end];
    NN=unique(NN);
    
    tic
    
    CDF_min =0;
    
    while CDF_min <= epsil %continue until we find the correct eta such that CDF_min>epsil
        
        total_samp_num = length(NN)*num_loop_x1*(num_loop_x3-1);
        
        CDF_array=-1*ones( num_loop_x1*(num_loop_x3-1),length(NN));
        CWD_array=zeros( num_loop_x1*(num_loop_x3-1)*3, length(NN));
        
        parfor index_N = 1: length(NN)
            
            N=NN(index_N);
            x1_min  = E/(1+N); %in this case x1= x2, x3=0, this is to make sure that x1>=x2>=x3.    
            step_x1 = E/num_loop_x1; 
            x1_array = [x1_min : step_x1 : (E-step_x1)]; % the length of x1_array is at most num_loop_x1
            CDF_temp=-1*ones(num_loop_x1*(num_loop_x3-1),1);
            CWD_temp=zeros(3,num_loop_x1*(num_loop_x3-1));
            
            for index_x1 = 1: length(x1_array)
                
                x1 = x1_array(index_x1); 
                x3_max = (E-x1)/(N+1); % the max of x3 is equal to x2
                x3_array=[1:num_loop_x3-1]*x3_max/num_loop_x3;% the length of x3_array is num_loop_x3-1 
               
                
                for index_x3 = 1:length(x3_array) 
                    
                    x3=x3_array(index_x3);
                    
                    x2= (E-x1-x3)/N;
                    
                    if x2>=x1-eps_small | x3>=x2 -eps_small | x3<= eps_small
                        continue
                    end
                    
                    R_prob = eta + log(1+x1) + N*log(1+x2) + log(1+x3);
                    
                    CDF = gammainc(R_prob/x2,N) - x1/(x1-x3)* exp(-R_prob/x1 - N*log(1 -x2/x1) )* gammainc(R_prob*(1/x2-1/x1),N);
                    arg_int=0: R_prob/10000:R_prob;
                    arg_int(1)=arg_int(1)+1e-16;
                    pdf_GN = exp((N-1)*log(arg_int) - arg_int/x2 - N*log(x2) - gammaln(N) );
                    
                    int2= sum( exp(-(R_prob-arg_int)/x3).* pdf_GN  )*(R_prob/10000);
                    
                    CDF = CDF + x3/(x1-x3)* int2;
                    
                    CDF_temp((index_x1-1)*(num_loop_x3-1) + index_x3) = CDF;
                    
                    CWD_temp(:,(index_x1-1)*(num_loop_x3-1) + index_x3)= [x1,x3,N].';
                end % end x3
            end %end x1
            
            CDF_array(:,index_N)=CDF_temp;
            CWD_array(:,index_N)= reshape(CWD_temp,3*num_loop_x1*(num_loop_x3-1),1);
            
        end %end N
        
        
        
        %Next, we consider the case where the optimal codeword contains two
        %distinct nonzero values 
        step3 = (log(N_end) -log(N_start))/num_loop_NN2;
        NN2 = [round(exp(log(N_start):step3: log(N_end))), N_end];
        NN2=unique(NN2);
        
        total_samp_num2 = length(NN2)^2*num_loop_x12;
        CDF_array2=-1*ones( length(NN2)*num_loop_x12, length(NN2));
        CWD_array2=zeros(3*length(NN2)*num_loop_x12, length(NN2));
        
        parfor index_N1=1:length(NN2)
            
            N1 = NN2(index_N1);
            CDF_temp = -1*ones(length(NN2)*num_loop_x12,1);
            CWD_temp = zeros(3,length(NN2)*num_loop_x12);
            
            for index_N2=1:length(NN2)
                
                N2 =NN2(index_N2);    
                step_x1= (E/N1 -E/(N1+N2))/(num_loop_x12-1);
                x1_array = E/(N1+N2): step_x1: E/N1;
                
                for index_x1 = 1:length(x1_array)
                    
                    x1 = x1_array(index_x1);
                    x2 = (E-x1*N1)/N2;
                    if x2<eps_small | x2>=x1-eps_small
                        continue;
                    end
                    R_prob = eta + N1*log(1+x1) + N2*log(1+x2);
                    arg_int=0: R_prob/10000:R_prob;
                    PDF_G2 = exp((N2-1)*log(arg_int) - arg_int/x2 - N2*log(x2) - gammaln(N2));
                    CDF_G1 = gammainc((R_prob-arg_int)/x1, N1);
                    CDF = sum(PDF_G2.*CDF_G1)*(R_prob/10000);
                    
                    CDF_temp((index_N2-1)*num_loop_x12+index_x1,1)=CDF;
                    CWD_temp(:,(index_N2-1)*num_loop_x12+index_x1)=[x1,N1,N2].';
                    
                end
            end
            CDF_array2(:,index_N1) = CDF_temp;
            CWD_array2(:,index_N1) = reshape(CWD_temp, 3*num_loop_x12*length(NN2),1);
        end
        
        CDF_array=reshape(CDF_array, 1, num_loop_x1*(num_loop_x3-1)*length(NN));
        CWD_array=reshape(CWD_array, 3, num_loop_x1*(num_loop_x3-1)*length(NN)).';
        
        indexx = find(CDF_array<=0);
        CDF_array(indexx)=[];
        CWD_array(indexx,:)=[];
        
        CDF_array2=reshape(CDF_array2, 1, length(NN2)^2*num_loop_x12);
        CWD_array2=reshape(CWD_array2, 3, length(NN2)^2*num_loop_x12).';
        
        indexx = find(CDF_array2<=0);
        CDF_array2(indexx)=[];
        CWD_array2(indexx,:)=[];
        
        
        CDF_min = min([min(CDF_array), target_error,min(CDF_array2)]);
        
        %We need to update eta if CDF_min<= epsil
        
        if CDF_min<=epsil  
            
            
            
            if min(CDF_array) < min(CDF_array2) % the optimal codeword contains three distinct nonzero values
                
                opt_cwd_index = find(CDF_array == min(CDF_array));
                x1_opt = CWD_array(opt_cwd_index(1),1);
                x3_opt = CWD_array(opt_cwd_index(1),2);
                N_opt  = CWD_array(opt_cwd_index(1),3);
                x2_opt= (E-x1_opt-x3_opt)/N_opt;
                
                R_prob = eta + log(1+x1_opt) + N_opt*log(1+x2_opt) + log(1+x3_opt);
                CDF =0;
                while CDF <= target_error
                    
                    R_prob = R_prob*1.05;
                    
                    CDF = gammainc(R_prob/x2_opt,N_opt) - x1_opt/(x1_opt-x3_opt)* exp(-R_prob/x1_opt - N_opt*log(1 -x2_opt/x1_opt) )* gammainc(R_prob*(1/x2_opt-1/x1_opt),N_opt);
                    
                    arg_int=0: R_prob/10000:R_prob;
                    
                    arg_int(1)=arg_int(1)+1e-16;
                    
                    pdf_GN = exp((N_opt-1)*log(arg_int) - arg_int/x2_opt - N_opt*log(x2_opt) - gammaln(N_opt) );
                    
                    int2= sum( exp(-(R_prob-arg_int)/x3_opt).* pdf_GN  )*(R_prob/10000);
                    
                    CDF = CDF + x3_opt/(x1_opt-x3_opt)* int2;
                end
                
                eta = R_prob -log(1+x1_opt) - N_opt*log(1+x2_opt) - log(1+x3_opt);
                
            else% the optimal codeword contains two distinct nonzero values
                
                opt_cwd_index = find(CDF_array2 ==min(CDF_array2));
                x1_opt = CWD_array2(opt_cwd_index(1),1);
                N1_opt = CWD_array2(opt_cwd_index(1),2);
                N2_opt = CWD_array2(opt_cwd_index(1),3);
                x2_opt = (E-x1_opt*N1_opt)/N2_opt;
                
                R_prob = eta + N1_opt*log(1+x1_opt) + N2_opt*log(1+x2_opt);
                CDF=0;
                
                while CDF<=target_error
                    R_prob = R_prob*1.05;
                    arg_int=0: R_prob/10000:R_prob;
                    PDF_G2 = exp((N2_opt-1)*log(arg_int) - arg_int/x2_opt - N2_opt*log(x2_opt) - gammaln(N2_opt));
                    CDF_G1 = gammainc((R_prob-arg_int)/x1_opt, N1_opt);
                    
                    CDF = sum(PDF_G2.*CDF_G1)*(R_prob/10000);
                end
                eta = R_prob - N1_opt*log(1+x1_opt) - N2_opt*log(1+x2_opt); 
            end   
        end %end if CDF_min<=epsil
            
    end%end while CDF_min<=epsil
     
    info_k =  (eta - log(CDF_min -epsil))/log(2);
    info_bits = [info_bits, info_k];
    toc
    
    
end