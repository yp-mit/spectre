function [R,current_eps]=MC_USTM_Alamouti(snrdB,T,L,epsilon,prec,filename)
%
% metaconverse bound for the Alamouti ensemble
%-------------------------------------------------------------------
%                       SET-UP PARAMETERS
%-------------------------------------------------------------------

SAVE=1;
MAT=1;

K = 2^prec; % number of monte carlo simulations (at least 100 x 1/epsilon)
rho = 10.^(snrdB/10); % SNR in linear scale
Mt=2;
Mr=2;
%-------------------------------------------------------------------
%                       MONTE CARLO SIMULATION
%-------------------------------------------------------------------
Ip = zeros(K,1); %allocate for the montecarlo runs
% Iq = zeros(K,1); %allocate for the montecarlo runs
%-------------------------------------------------------------------
%                       CONSTANTS
%-------------------------------------------------------------------
rho_tilde = T*rho/Mt;

lambda=1+rho_tilde;
lambda1=1/lambda;
lambda2=rho_tilde*lambda1;

x1=rho_tilde;
x2=rho_tilde;

D= [diag([sqrt(1+x1), sqrt(1+x2)]), zeros(Mt,T-Mt);
zeros(T-Mt,Mt), eye(T-Mt)]; % D matrix (covariance matrix of equivalent noise)


%-------------------------------------------------------------------
%                       MONTE CARLO
%-------------------------------------------------------------------

%tic


norm=sqrt(.5);

for k=1:K
    

        i_L = 0; 
        
        Z = randn(T,Mr,L)*norm+1i*randn(T,Mr,L)*norm;  

        for l = 1:L %Create each realization
                  
            
          %COMPUTE EVERYTHING THAT HAS TO DO WITH SINGULAR VALUES
          Sigma_alt= svd(D*Z(:,:,l)).^2;
          Sigma_alt=sort(Sigma_alt,1,'descend');
          TraceZ=abs(trace(Z(:,:,l)'*Z(:,:,l)));
          
          Y=D*Z(:,:,l);
          
          Ytilde=(zeros(T,4)); 
          
          Ytilde(:,[1,3])=Y;
          
          Ytilde(1:2:T,[2,4])= conj(Y(2:2:T,:));
          
          Ytilde(2:2:T,[2,4])= -conj(Y(1:2:T,:));
          
          Sigma=svd(Ytilde).^2;
          
          Sigma=[Sigma(1),Sigma(3)];
          
          Sigma=Sigma*lambda2;
      
        
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
   
        end

        Ip(k) =  i_L; %put all computations on a pile to compute the average later
end

if (SAVE==1) 
  if (MAT==1)
    save(filename,'Ip')
  else
    save(filename,'Ip','-ascii','-append')
  end
end



%---------------------------------------
%   SEARCHING THE RATE
%--------------------------------------- 

% load saved data (to account for append possibilities)
if (SAVE==1 && MAT==0)

    Ip=load(filename);
end

Ip=sort(Ip);

Kcurrent=length(Ip); % redefine K to account for append

current_prec=floor(log2(Kcurrent)); % actual precision

K=2^(current_prec); % round off K to avoid search errors


% first find a suitable initial point for the linear search
step=K/2;
index=step;

onevec=ones(K,1);

while(step>1),
    
   th=Ip(index);
   
   current_eps=sum(Ip<=th)/K;
   
   step=step/2;
   
   if current_eps> epsilon,
       
       index=index-step;
       
   else
       
       index=index+step;
       
   end
   
end

current_eps=sum(Ip<=Ip(index))/K; 

if(current_eps<epsilon)
    index=index+1;
end





% now perform linear search
% Rvect=zeros(1,K-index+1);
% for ii=1:K-index+1,
% %for ii=1:K-index+1,
%
%     current_rate=Ip(ii+index-1)-log(sum(Ip<=Ip(ii+index-1))/K-epsilon);
%
%     Rvect(ii)=current_rate;
%
% end

%faster search (potentially less tight)

count=0;

R=Ip(index)-log(sum(Ip<=Ip(index))/K-epsilon);

for ii=1:K-index+1, 
    current_rate=Ip(ii+index-1)-log(sum(Ip<=Ip(ii+index-1))/K-epsilon);  
    if current_rate<= R
        R=current_rate;
     else
       count=count+1;
    end
     
    if count==20,
      break
    end
end
R=R/(L*T*log(2)); 

    

end


