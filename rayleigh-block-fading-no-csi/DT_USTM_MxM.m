
function [R,current_eps]=DT_USTM_MxM(snrdB,T,L,Mant,epsilon,prec,filename)
%
% Function to compute the DT lower bound for a Rayleigh block-fading
% channel with no CSI at transmitter and receiver. The bound assumes that
% USTM is chosen as input distribution 
% snrdB: SNR in dB
% Mant: number of transmit antennas and number of receive antennas (square setting)
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
CLUSTER=0;
MAT=0;

Mt=Mant; %number of transmit antennas
Mr=Mant; %number of receive antennas
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
% D= [(1+rho_tilde)*eye(Mt), zeros(Mt,T-Mt);
%     zeros(T-Mt,Mt), eye(T-Mt)]; % D matrix (covariance matrix of equivalent noise)
D= [sqrt(1+rho_tilde)*eye(Mt), zeros(Mt,T-Mt);
    zeros(T-Mt,Mt), eye(T-Mt)]; % D matrix (covariance matrix of equivalent noise)



lambda=1+rho_tilde;
lambda1=1/lambda;
lambda2=rho_tilde*lambda1;
c2 = logComplexGammaRatio(Mt,Mr,T); %gamma constant
c1 = Mt*(T-Mt)*log(lambda2); % SNR constant
const = c1+c2;

P=max(Mt,Mr);


%-------------------------------------------------------------------
%                       MONTE CARLO
%-------------------------------------------------------------------

if (CLUSTER==1),
    [~, tmpdir] = system('echo $TMPDIR');
    matlabpool close force;
    sched = findResource('scheduler', 'configuration', 'local');
    sched.DataLocation = tmpdir(1:end-1);
    matlabpool open local;
end  

tic

for k = 1:K %do K montecarlo runs 

        i_L = 0;   
  
        for l = 1:L %Create each realization
            %GENERATE Z
            Z = randn(T,Mr)*sqrt(.5)+1i*randn(T,Mr)*sqrt(.5);
            Sigma=svd(D*Z).^2; 
            M = createM(Mt,Mr,P,T,Sigma,lambda2); %create  matrix
            detM = det(M); %get the determinant
            vanderTerm = det(vander(Sigma)); %get the determinant of the vandermode matrix
            TraceZ=trace(Z'*Z);      
            logDetSigma = (T-Mr)*sum(log(Sigma));% compute det(Sigma^(T-M))
            logdetM=log(detM);
            i = const- TraceZ  + lambda1*sum(Sigma) - log(detM) + log(vanderTerm) + logDetSigma;
            i_L = i_L + i; %add it to the total i_L 
            
        end

        I(k) =  i_L; %put all computations on a pile to compute the average later
        
end

actual_time=toc/60;

msg=sprintf('time needed to generate samples: %f min', actual_time);
disp(msg); 

%disp('samples generated')

if (SAVE==1) 
  if (MAT==1)
    save(filename,'I')
  else
    save(filename,'I','-ascii','-append')
  end
end

%disp('samples saved')

if (CLUSTER==1),
    matlabpool close;
end

%---------------------------------------
%   START SEARCHING FOR THE RATE
%---------------------------------------  
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


R=log2(exp(I(index))+1)/(L*T); % factor 2 removed from DT bound to account for max error probability

end


function val = logComplexGammaRatio(Mt,Mr,T)
    k=T-min(Mt,Mr)+1:1:T;
    val = -sum(gammaln(k));
end

function M = createM(Mt,Mr,P,T,Sigma,lambda)
    M = nan(P,P); %(l is the row, k is the column)
    % note: this function is invoked in the following program only for the case Mt=Mr;
    for l = 1:Mr,
        for k=1:Mt,
                
            M(l,k)=(Sigma(l)^(Mt-k))*gammainc(lambda*Sigma(l),T+k-Mt-P);
        end
        
    end

                
end

 





