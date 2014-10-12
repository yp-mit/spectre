function [R,current_eps,current_prec]=MC_USTM_2x2_Qmimo(snrdB,T,L,M,epsilon,prec,pow_all,filename)
%
% Metaconverse upper bound based on an auxiliary output distribution induced by a 2xT USTM input distribution
%
%
% %
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
%   lpow_all: L-dimensional power allocation vector: 0.5 means equal power on
% each antenna; 0 means no power on the second antenna; Note that each entry of pow_all can be taken
% to be between 0 and 0.5 without loss of generality.
%
%   filename: data file where the samples of the information density are saved for possible future refinements. 
%
% The outputs of the program are
% R: R^*(T,L,current_eps,snrdB)
% current_eps: actual upper bound on the error probability for the given
% code
% current_prec: returns log2 of the number of samples used to determine R.
%-------------------------------------------------------------------
%                       SET-UP PARAMETERS
%-------------------------------------------------------------------

SAVE=1; % append samples of information density to the file "filename"


K = 2^prec; % number of monte carlo simulations (at least 100 x 1/epsilon)
rho = 10.^(snrdB/10); % SNR in linear scale
Mt=M;
Mr=M;
%-------------------------------------------------------------------
%                       MONTE CARLO SIMULATION
%-------------------------------------------------------------------
Ip = zeros(K,1); %allocate for the montecarlo runs
%-------------------------------------------------------------------
%                       CONSTANTS
%-------------------------------------------------------------------
rho_tilde = T*rho/Mt;


lambda=1+rho_tilde;
lambda1=1/lambda;
lambda2=rho_tilde*lambda1;
c2 = logComplexGammaRatio(Mt,Mr,T); %gamma constant

P=max(Mt,Mr);

%-------------------------------------------------------------------
%                       MONTE CARLO
%-------------------------------------------------------------------


for k=1:K
    

        i_Lp = 0;   

        for l = 1:L %Create each realization
            
            x1=rho_tilde*(2*pow_all(l));
            x2=rho_tilde*2*(1-pow_all(l));

            D= [diag([1+x1, 1+x2]), zeros(Mt,T-Mt);
    zeros(T-Mt,Mt), eye(T-Mt)]; % D matrix (covariance matrix of equivalent noise)

            
            c1 = Mt*(T-Mt)*log(lambda2)+Mr*log(lambda/(1+x1))+Mr*log(lambda/(1+x2)); % SNR constant
            const = c1+c2;

            
            % compute samples of information density log dPy|x/dQy for y~ Py|x
            Z = randn(T,Mr)*sqrt(.5)+1i*randn(T,Mr)*sqrt(.5);
            Sigma = abs(eig(Z'*D*Z));
            Sigma=sort(Sigma,1,'descend');
            M = createM(Mt,Mr,P,T,Sigma,lambda2); %create  matrix
            detM = det(M); %get the determinant
            vanderTerm = det(vander(Sigma)); %get the determinant of the vandermode matrix
            logdetM=log(detM);

            TraceZ=trace(Z'*Z);
            logDetSigma = (T-Mr)*sum(log(Sigma));% compute det(Sigma^(T-M))
            
            ip = const- TraceZ  + lambda1*sum(Sigma) - logdetM + log(vanderTerm) + logDetSigma;%Information density for time l (i(x_l, y_l)
            i_Lp = i_Lp + ip; %add it to the total i_L 
                        
   
        end

        Ip(k) =  i_Lp; %put all computations on a pile to compute the average later
end



if (SAVE==1)

    save(filename,'Ip','-ascii','-append')
end


%---------------------------------------
%   SEARCHING THE RATE
%--------------------------------------- 

% load saved data (to account for append possibilities)
if (SAVE==1)

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
Rvect=zeros(1,K-index+1);


for ii=1:K-index+1,
    
    current_rate=Ip(ii+index-1)-log(sum(Ip<=Ip(ii+index-1))/K-epsilon);
    
    Rvect(ii)=current_rate;
    
end
    


R=min(Rvect)/(L*T*log(2));

end


function val = logComplexGammaRatio(Mt,Mr,T)
    k=T-min(Mt,Mr)+1:1:T;
    val = -sum(gammaln(k));
end

function M = createM(Mt,Mr,P,T,Sigma,lambda)
    M = nan(P,P); %(l is the row, k is the column)

    
    for l = 1:Mr,
        for k=1:Mt,
                
            M(l,k)=(Sigma(l)^(Mt-k))*gammainc(lambda*Sigma(l),T+k-Mt-P);
        end
        
        for k=Mt+1:Mr,
                    
            M(l,k)=(Sigma(l)^(T-k))*exp(-Sigma(l)*lambda); % case Mt<Mr
            
        end
    end

                
end