function [Cout,current_eps]=outage_capacity_alamouti(snrdb,epsilon,tfdiv,prec,filename)
%
%
% compute outage capacity of a 2x2 block-memoryless MIMO Rayleigh fading
% channel when an Alamouti inner code is used
% snrdb: snr in dB
% epsilon: outage probability
% tfdiv: number of time_frequency diversity branches
% prec: it controls the number of samples for the Monte Carlo simulation; Note nsamples=2^prec; One should have nsamples>> 100 x 1/epsilon


%%%%%%%%%%%%%%%%%%%%%%%%
% SET-UP parameters
%%%%%%%%%%%%%%%%%%%%%%%%

SAVE=1;

mt=2; %number of transmit antennas

mr=2; %number of receive antennas

K = 2^prec; % number of monte carlo simulations (it should be at least 100 x 1/epsilon)
rho = (10.^(snrdb/10))/mt; % SNR in linear scale



norm=sqrt(0.5);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONTECARLO: generation of samples of instantaneous mutual information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



mi=zeros(K,1);

for ii=1:K,
    
    H=(randn(mt,mr,tfdiv)+1i*randn(mt,mr,tfdiv))*norm;
    
    for k=1:tfdiv
    
        mi(ii)=mi(ii)+log2(1+rho*abs(trace(H(:,:,k)*H(:,:,k)')));
       
    end
    
end

mi=mi/tfdiv;



if (SAVE==1)
    save(filename,'mi','-ascii','-append')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation outage capacity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (SAVE==1)
    mi=load(filename);
end

mi=sort(mi,'ascend');

step=K/2;
index=step;

onevec=ones(K,1);

while(step>1),
    
   th=mi(index);
   
   current_eps=sum(mi<=th)/K;
   
   step=step/2;
   
   if current_eps> epsilon,
       
       index=index-step;
       
   else
       
       index=index+step;
       
   end
   
end

current_eps=sum(mi<=mi(index))/K;

if(current_eps<epsilon)
    index=index+1;
end
current_eps=sum(mi<=mi(index))/K;

Cout=mi(index);