
function [R,current_eps]=DT_USTM_Alamouti(snrdB,T,L,epsilon,prec,filename)
%
% Function to compute the DT lower bound for a 2x2 Rayleigh block-fading
% channel with no CSI at transmitter and receiver. The bound assumes that
% USTM+Alamouti innner code is chosen as input distribution
% snrdB: SNR in dB
% T: size of coherence interval
% L: number of independent coherence intervals
% epsilon: average block error probability
%
% The outputs of the program are
% R: R^*(T,L,current_eps,snrdB)
% current_eps: actual upper bound on the error probability for the given
% code
%-------------------------------------------------------------------
%                       SET-UP PARAMETERS
%-------------------------------------------------------------------

SAVE=1;
MAT=1;
CLUSTER=1;


K=2^prec;

rho = 10.^(snrdB/10); % SNR in linear scale

%     N = M; % receive antennas 
%     Lmax = n/(N+M); %find the largest number of coherence times to simulate
%     L =  1:Lmax; %all # coherence times of interest
%     Ttmp = n./L; %get corresponding coherence times
%     indx = find((rem(Ttmp,1) == 0) == 1); %only take integer T
%     L=L(indx);
%     coherence_times = n./L;% coherence times
%-------------------------------------------------------------------
%                       MONTE CARLO SIMULATION
%-------------------------------------------------------------------
rates = nan(1); %allocate vector to save rates
I = zeros(K,1); %allocate for the montecarlo runs
%-------------------------------------------------------------------
%                       CONSTANTS
%-------------------------------------------------------------------
Mt=2;
Mr=2;
rho_tilde = T*rho/Mt; 
D= [sqrt(1+rho_tilde)*eye(Mt), zeros(Mt,T-Mt);
    zeros(T-Mt,Mt), eye(T-Mt)]; % D matrix (covariance matrix of equivalent noise)

lambda=1+rho_tilde;
lambda1=1/lambda;
lambda2=rho_tilde*lambda1; % beta in my notes

norm=sqrt(.5);

%-------------------------------------------------------------------
%                       MONTE CARLO
%-------------------------------------------------------------------

tic

if (CLUSTER==1),
    [~, tmpdir] = system('echo $TMPDIR');
    matlabpool close force;
    sched = findResource('scheduler', 'configuration', 'local');
    sched.DataLocation = tmpdir(1:end-1);
    matlabpool open local;
end
    
parfor k = 1:K %do K montecarlo runs

        Z = randn(T,Mr,L)*norm+1i*randn(T,Mr,L)*norm;
        
        i_L = 0;   
        for l = 1:L %Create each realization
            %GENERATE Z
            
            %COMPUTE EVERYTHING THAT HAS TO DO WITH SINGULAR VALUES
            Sigma_alt= svd(D*Z(:,:,l)).^2;
            TraceZ=abs(trace(Z(:,:,l)'*Z(:,:,l)));
            
            Y=D*Z(:,:,l);
            
            Ytilde=(zeros(T,4)); 
            
            Ytilde(:,[1,3])=Y;
            
            Ytilde(1:2:T,[2,4])= conj(Y(2:2:T,:));
            
            Ytilde(2:2:T,[2,4])= -conj(Y(1:2:T,:));
            
            Sigma=svd(Ytilde).^2;
            
            Sigma=[Sigma(1),Sigma(3)];
            
            Sigma=Sigma*lambda2;
            
            % now compute the singular values of the Y_tilde matrix and keep only first and third
            
            %%% Wei's implementation
            % M1=[gammainc(Sigma(1), T-5,'scaledlower')/Sigma(1) ,  gammainc(Sigma(1), T-4,'scaledlower')/(T-4), gammainc(Sigma(2), T-5,'scaledlower')/Sigma(2), gammainc(Sigma(2), T-4,'scaledlower')/(T-4) ;
   %                          (T-2)*Sigma(1), Sigma(1)^2,(T-2)*Sigma(2), Sigma(2)^2 ;
   %                          T-3, Sigma(1),T-3, Sigma(2);
   %                          (T-4)/Sigma(1),1,(T-4)/Sigma(2),1];
   %
   %          d1=det(M1);
   %
   %          log_exp_sum=log(d1)- 4*log(Sigma(1)-Sigma(2));  
   
   %          i = - TraceZ  +sum(Sigma_alt) - log(T-1)-log(T-2)-log(T-3)-log(T-4) - log_exp_sum; 
               

            %%% Giuseppe's implementation 1
            
            % norm=1;
     %
     %        M=[exp(Sigma(1)-norm*Sigma(1))*gammainc(Sigma(1), T-5)/Sigma(1)^(T-4) ,  exp(Sigma(1)-norm*Sigma(1))*gammainc(Sigma(1), T-4)/Sigma(1)^(T-4), exp(Sigma(2)-norm*Sigma(1))*gammainc(Sigma(2), T-5)/Sigma(2)^(T-4), exp(Sigma(2)-norm*Sigma(1))*gammainc(Sigma(2), T-4)/Sigma(2)^(T-4) ;
     %                (T-2)*Sigma(1), Sigma(1)^2,(T-2)*Sigma(2), Sigma(2)^2 ;
     %                T-3, Sigma(1),T-3, Sigma(2);
     %                (T-4)/Sigma(1),1,(T-4)/Sigma(2),1]
     
             %log_exp_sum = log(d) -4*log(Sigma(1)-Sigma(2))-(T-4)*log(Sigma(1))+Sigma(1);
   
             % Giuseppe's implementation 2
             
          if (T>4),  
               
            M=[gammainc(Sigma(1), T-5),  gammainc(Sigma(1), T-4), exp(Sigma(2)-Sigma(1))*gammainc(Sigma(2), T-5)/(Sigma(2)/Sigma(1))^(T-4), exp(Sigma(2)-Sigma(1))*gammainc(Sigma(2), T-4)/(Sigma(2)/Sigma(1))^(T-4) ;
             (T-2)*Sigma(1), Sigma(1)^2,(T-2)*Sigma(2), Sigma(2)^2 ;
             T-3, Sigma(1),T-3, Sigma(2);
             (T-4)/Sigma(1),1,(T-4)/Sigma(2),1];
            
            logd=log(det(M))-(T-4)*log(Sigma(1))+Sigma(1);
            
          else
            
            M=[ 1, 2*Sigma(1), 1, 0 ; ...
                1,  Sigma(1)^2, Sigma(1), 1; ...
                exp(Sigma(2)-Sigma(1)), 2*Sigma(2), 1, 0; ...
                exp(Sigma(2)-Sigma(1)), Sigma(2)^2, Sigma(2), 1];
             
            logd=(log(det(M)))+Sigma(1);
                
          end
              
            
            log_exp_sum = logd -4*log(Sigma(1)-Sigma(2));
          
          
          
            i = - TraceZ  +sum(Sigma_alt) - log(gamma(T)) - log_exp_sum;
            
            i_L = i_L + i; %add it to the total i_L 
            
%             if isnan(i)
%                 sprintf('i error')
%             end
%             if isnan(det(M))
%                 sprintf('error')
%             end
        end

        I(k) =  i_L; %put all computations on a pile to compute the average later
        
        
end

actual_time=toc/60;

msg=sprintf('time needed to generate samples: %f min', actual_time);
disp(msg); 

if (CLUSTER==0)
  est_time=actual_time*2^(23-prec)/6;
  msg=sprintf('expected cluster time: %f min', est_time);
  disp(msg);
end
  


if (SAVE==1) 
  if (MAT==1)
    save(filename,'I')
  else
    save(filename,'I','-ascii','-append')
  end
end

if (CLUSTER==1),
    matlabpool close;
end


%---------------------------------------
%   START SEARCHING FOR THE RATE
%---------------------------------------  
% load saved data (to account for append possibilities)
if (SAVE==1 && MAT==0)

    I=load(filename);
end

I=sort(I);

Kcurrent=length(I); % redefine K to account for append

current_prec=floor(log2(Kcurrent)); % actual precision

K=2^(current_prec); % round off K to avoid search errors



step=K/2;
index=step;

onevec=ones(K,1);

while(step>1),
    
   th=I(index);
   
   current_eps=sum(exp(-max(0,I-th)))/K;
   
   step=step/2;
   
   if current_eps> epsilon,
       
       index=index-step;
       
   else
       
       index=index+step;
       
   end
   
end

current_eps=sum(exp(-max(0,I-I(index))))/K;


R=log2(2*exp(I(index))+1)/(L*T);

end




 





