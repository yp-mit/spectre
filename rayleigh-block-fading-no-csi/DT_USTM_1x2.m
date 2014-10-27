function [R,current_eps]=DT_USTM_1x2(snrdB,T,L,epsilon,prec,filename)
% DT_USTM_1x2 description
%   Function to compute the DT lower bound for a Rayleigh block-fading
%   channel with no CSI at transmitter and receiver. The bound assumes that
%   USTM is chosen as input distribution, that 1 transmit antenna is used and two receive antennas are available at the receiver
%
%   snrdB: SNR in dB
%
%   T: size of coherence interval
%
%   L: number of independent coherence intervals
%
%   epsilon: average block error probability
%
%   prec: it controls the number of samples for the Monte Carlo simulation; Note nsamples=2^prec; One should have nsamples>> 100 x 1/epsilon
%
%   filename: data file where the samples of the information density are saved for possible future refinements. 
%
% The outputs of the program are
% R: R^*(T,L,current_eps,snrdB)
% current_eps: actual upper bound on the error probability for the given
% code
%-------------------------------------------------------------------
%                       SET-UP PARAMETERS
%-------------------------------------------------------------------

SAVE=1;

K=2^prec;

rho = 10.^(snrdB/10); % SNR in linear scale

%-------------------------------------------------------------------
%                       MONTE CARLO SIMULATION
%-------------------------------------------------------------------
rates = nan(1); %allocate vector to save rates
I = zeros(K,1); %allocate for the montecarlo runs
%-------------------------------------------------------------------
%                       CONSTANTS
%-------------------------------------------------------------------
Mt=1;
Mr=2;
rho_tilde = T*rho/Mt; 
D= [(1+rho_tilde)*eye(Mt), zeros(Mt,T-Mt);
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

    
for k = 1:K %do K montecarlo runs
    
    %waitbar(k/K)

        i_L = 0;   
        for l = 1:L %Create each realization
            %GENERATE Z
            Z = randn(T,Mr)*sqrt(.5)+1i*randn(T,Mr)*sqrt(.5);
            
            %COMPUTE EVERYTHING THAT HAS TO DO WITH SINGULAR VALUES
            Sigma = abs(eig(Z'*D*Z));
            Sigma=sort(Sigma,1,'descend');
            vanderTerm = det(vander(Sigma)); %get the determinant of the vandermode matrix
            TraceZ=trace(Z'*Z);
            
            a1=gammainc(lambda2*Sigma(1),T-2);
            a2=gammainc(lambda2*Sigma(2),T-2);
            
            sigmaprod=(Sigma(1)/Sigma(2))^(T-2);
            
            expratio=exp(-lambda2*(Sigma(1)-Sigma(2)));
            
            a3= (a2/a1)*sigmaprod * expratio;
         
            logdetM=log(a1) +(T-2)*log(Sigma(2))-Sigma(2)*lambda2+log( 1- a3);
        
            logDetSigma = (T-Mr)*sum(log(Sigma));% compute det(Sigma^(T-M))

            % GET INFORMATION DENSITY
            i = const- TraceZ  + lambda1*sum(Sigma) - logdetM + log(vanderTerm) + logDetSigma;%Information density for time l (i(x_l, y_l)
            i_L = i_L + i; %add it to the total i_L 
         
        I(k) =  i_L; %put all computations on a pile to compute the error probability later
end



if (SAVE==1)

    save(filename,'I','-ascii','-append') %save information density samples for future refinements
end


%---------------------------------------
%   START SEARCHING FOR THE RATE
%---------------------------------------  
if (SAVE==1)

    I=load(filename); % load information density samples
end

I=sort(I);

Kcurrent=length(I); % redefine K to account for append operation

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


function val = logComplexGammaRatio(Mt,Mr,T)
    k=T-min(Mt,Mr)+1:1:T;
    val = -sum(gammaln(k));
end



 





