function [R,current_eps,current_prec]=MC_USTM_2x2_q2x2(snrdB,T,L,Mt,epsilon,prec,pow_all,filename)
%

%-------------------------------------------------------------------
%                       SET-UP PARAMETERS
%-------------------------------------------------------------------
SAVE=1;
MAT=0; % if set, the .mat extension is used to save the file

K = 2^prec; % number of monte carlo simulations (at least 100 x 1/epsilon)
rho = 10.^(snrdB/10); % SNR in linear scale
Mr=Mt;
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
c2 = logComplexGammaRatio(Mt,Mr,T); %gamma constant

P=max(Mt,Mr);

%-------------------------------------------------------------------
%                       MONTE CARLO
%-------------------------------------------------------------------

noise_norm=sqrt(.5);

%parfor k = 1:K %do K montecarlo runs
% tic
parfor k=1:K
    

        i_Lp = 0;  
        
        Z = randn(T,Mr,L)*noise_norm+1i*randn(T,Mr,L)*noise_norm; 

        for l = 1:L %Create each realization
            
            x1=rho_tilde*(2*pow_all(l));
            x2=rho_tilde*2*(1-pow_all(l));

            D= [diag([sqrt(1+x1), sqrt(1+x2)]), zeros(Mt,T-Mt);
    zeros(T-Mt,Mt), eye(T-Mt)]; % D matrix (covariance matrix of equivalent noise)
            
            c1 = Mt*(T-Mt)*log(lambda2)+Mr*log(lambda/(1+x1))+Mr*log(lambda/(1+x2)); % SNR constant
            const = c1+c2;

            
            % compute samples of information density log dPy|x/dQy for y~ Py|x
            
            Sigma=svd(D*Z(:,:,l)).^2; 
            SigmaAlt=Sigma*lambda2;
           % M = createM(Mt,Mr,P,T,Sigma,lambda2); %create  matrix
            
            M=[Sigma(1)*gammainc(SigmaAlt(1),T-3), gammainc(SigmaAlt(1),T-2); ...
            Sigma(2)*gammainc(SigmaAlt(2),T-3), gammainc(SigmaAlt(2),T-2)];
            
            detM = det(M); %get the determinant
            vanderTerm = det(vander(Sigma)); %get the determinant of the vandermode matrix
            logdetM=log(detM);

            TraceZ=abs(trace(Z(:,:,l)'*Z(:,:,l)));
            logDetSigma = (T-Mr)*sum(log(Sigma));% compute det(Sigma^(T-M))
            
            ip = const- TraceZ  + lambda1*sum(Sigma) - logdetM + log(vanderTerm) + logDetSigma;%Information density for time l (i(x_l, y_l)
            i_Lp = i_Lp + ip; %add it to the total i_L 
                        
   
        end

        Ip(k) =  i_Lp; %put all computations on a pile to compute the average later
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
%Rvect=zeros(1,K-index+1);

% if (CLUSTER==1),
%     [~, tmpdir] = system('echo $TMPDIR');
%     matlabpool close force;
%     sched = findResource('scheduler', 'configuration', 'local');
%     sched.DataLocation = tmpdir(1:end-1);
%     matlabpool open local;
% end


%exact search
% for ii=1:K-index+1,
%
%     current_rate=Ip(ii+index-1)-log(sum(Ip<=Ip(ii+index-1))/K-epsilon);
%
%     Rvect(ii)=current_rate;
%
% end
% R=min(Rvect)/(L*T*log(2));
    
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