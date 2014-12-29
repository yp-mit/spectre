
function RergUSTM=ergodic_USTM(snrdB,T,Mt,Mr,prec)
%
% function to compute a lower bound on the ergodic capacity of a Rayleigh block-fading
% channel with no CSI at transmitter and receiver. The bound assumes that
% USTM is chosen as input distribution 
% snrdB: SNR in dB
% Mt: number of transmit antennas 
% Mr: number of receive antennas (Mt<=Mr)
% T: size of coherence interval 
% prec: it controls the number of samples for the Monte Carlo simulation; Note nsamples=2^prec; 
%
% The outputs of the program are
% RergUSTM: lower bound on the ergodic capacity with USTM inputs 
%
%-------------------------------------------------------------------
%                       SET-UP PARAMETERS
%-------------------------------------------------------------------


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
    

    Z = (randn(T,Mr)+1i*randn(T,Mr))*noise_norm;      
    Sigma = svd(D*Z).^2;  
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
    TraceZ=real(trace(Z'*Z));
    i_temp = const- TraceZ  +partial_sum + log(vanderTerm); %Information density 
    I(k) =  i_temp; %put all computations on a pile to compute the average later      
end

RergUSTM= mean(I)/(log(2)*T);



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

 





