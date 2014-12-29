
function [R,current_eps]=DT_USTM(snrdB,T,L,Mt,Mr,epsilon,prec,filename)
%
% Function to compute the DT lower bound for a Rayleigh block-fading
% channel with no CSI at transmitter and receiver. The bound assumes that
% USTM is chosen as input distribution 
% snrdB: SNR in dB
% Mt: number of transmit antennas 
% Mr: number of receive antennas (Mt<=Mr); 
% T: size of coherence interval 
% L: number of independent coherence intervals; L*T is the blocklength
% epsilon: maximal block error probability
% prec: it controls the number of samples for the Monte Carlo simulation; Note nsamples=2^prec; One should have nsamples>> 100 x 1/epsilon
% filename: data file where the samples of the information density are saved for possible future refinements. 
%
% The outputs of the program are
% R: R^*(T*L,current_eps,snrdB)
% current_eps: actual upper bound on the maximal error probability 
%
%-------------------------------------------------------------------
%                       SET-UP PARAMETERS
%-------------------------------------------------------------------

SAVE=0; % if this flag is active the samples of the information density are saved for possibe further refinemes
MAT=0; % save with .mat extension

K = 2^prec; % number of monte carlo simulations (it should be at least 100 x 1/epsilon)
rho = 10.^(snrdB/10); % SNR in linear scale

%-------------------------------------------------------------------
%                       MONTE CARLO SIMULATION
%-------------------------------------------------------------------
I = zeros(K,1); %allocate for the montecarlo runs
%-------------------------------------------------------------------
%                       CONSTANTS
%-------------------------------------------------------------------
rho_tilde = T*rho/Mt; 
D= [sqrt(1+rho_tilde)*eye(Mt), zeros(Mt,T-Mt);
    zeros(T-Mt,Mt), eye(T-Mt)]; % D matrix (covariance matrix of equivalent noise)

lambda=1+rho_tilde;
lambda1=1/lambda;
lambda2=rho_tilde*lambda1;
c2 = logComplexGammaRatio(Mt,Mr,T); %gamma constant
c1 = Mt*(T-Mt)*log(lambda2); % SNR constant
const = c1+c2;

noise_norm=sqrt(.5);


%-------------------------------------------------------------------
%                       MONTE CARLO
%-------------------------------------------------------------------

for k = 1:K %do K montecarlo runs
        i_L = 0;   
        Z = (randn(T,Mr,L)+1i*randn(T,Mr,L))*noise_norm;
        for l = 1:L %Create each realization
            Sigma = svd(D*Z(:,:,l)).^2;    
            if (Mt==1)
                M = createMalt(Mt,Mr,T,Sigma,lambda2); %create  matrix
                logdetM=log(det(M))+Sigma(1)*lambda2-(T-Mr)*log(Sigma(1)); 
                partial_sum=sum(Sigma)-logdetM;
            else
                M = createM(Mt,Mr,T,Sigma,lambda2);
                logdetM=log(det(M)); 
                logdetSigma = (T-Mr)*sum(log(Sigma));% compute logdet(Sigma^(T-M))
                partial_sum=lambda1*sum(Sigma) - logdetM +logdetSigma;
            end
            
            vanderTerm = det(vander(Sigma)); %get the determinant of the vandermode matrix
            TraceZ=real(trace(Z(:,:,l)'*Z(:,:,l)));
            
            i_temp = const- TraceZ  +partial_sum + log(vanderTerm);  %Information density for time l (i(x_l, y_l)
            i_L = i_L + i_temp; %add it to the total i_L 
        end

        I(k) =  i_L; %put all computations on a pile to compute the average later
        
end

 


if (SAVE==1) 
  if (MAT==1)
    save(filename,'I')
  else
    save(filename,'I','-ascii','-append')
  end
end



%---------------------------------------
%   START SEARCHING FOR THE RATE
%---------------------------------------  

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

if (I(index)>500),

  R=I(index)/(L*T*log(2));
  
else

  R=log2(exp(I(index))+1)/(L*T); % factor 2 removed from DT bound to account for max error probability

end

end


function val = logComplexGammaRatio(Mt,Mr,T)
    k=T-Mt+1:1:T;
    val = -sum(gammaln(k));
    r=1:1:Mt;
    val=val+sum(gammaln(r));
end

function M = createMalt(Mt,Mr,T,Sigma,lambda)
    M = nan(Mr,Mr); %(l is the row, k is the column)
    for l = 1:Mr,
      for k=1:Mt,      
        M(l,k)=(Sigma(l)^(Mt-k))*gammainc(lambda*Sigma(l),T+k-Mt-Mr)*exp(lambda*(Sigma(l)-Sigma(1)))*(Sigma(1)/Sigma(l))^(T-Mr);
      end
      
      for k=Mt+1:1:Mr,
          M(l,k)= Sigma(l)^(Mr-k);
      end
        
    end

                
end


function M = createM(Mt,Mr,T,Sigma,lambda)
    M = nan(Mr,Mr); %(l is the row, k is the column)
    P=Mr;
    
    for l = 1:Mr,
        for k=1:Mt,
                
            M(l,k)=(Sigma(l)^(Mt-k))*gammainc(lambda*Sigma(l),T+k-Mt-P);
        end
        
        for k=Mt+1:Mr,
                    
            M(l,k)=(Sigma(l)^(T-k))*exp(-Sigma(l)*lambda); % case Mt<Mr
            
        end
    end
    
end

 





