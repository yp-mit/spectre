%[rate_na , C_epsilon] = normapprox_mimo_iso(nn,P,error,tx,rx,K) computes
%the epsilon-capacity and the normal approximation for the maximal
%achievable rate over a quasi-static MIMO Rician fading channel with
%isotropic inputs; It also outputs the epsilon-capacity of the
%corresponding channel.    
%
%The function takes the following inputs:
%
%nn: blocklength, scalar/vector, please put nn in an ascending order
%P: total input power (linear scale), scalar;
%error: block error probability, scalar;
%tx: number of transmit antennas, scalar;
%rx: number of receive antennas, scalar;
%K:  rician K factor, scalar; the default value is 0.


function [rate_na, C_epsilon] = normapprox_mimo_iso(nn,P,error,tx,rx,K)
if (nargin < 6) || isempty(K)
	K = 0;
end
rate_na=[];
P=P/tx;

loop=round(1000/error);

C_array=zeros(1,loop);
V_array=zeros(1,loop);

min_tx_rx = min(tx,rx);
max_tx_rx = max(tx,rx);

loop1=10000;%precision of computation: 0.0001
rate_previous=100/loop1;


%Compute C(H,P) and V(H,P) for different realizations of H 
for ij=1:loop
    
    h = sqrt(0.5/(K+1))*(randn(min_tx_rx,max_tx_rx)+ 1i*randn(min_tx_rx,max_tx_rx))+sqrt(K/(K+1));    
    gain =abs((svd(h)).^2.)*P;
   
    %C_array(ij)=log(abs(det(eye(min_tx_rx) + P*h*h')));
    C_array(ij) = sum(log(1+gain));   
    V_array(ij) = min_tx_rx -sum (1./(1+gain).^2 );
end

C_s=sort(C_array);
%compute the epsilon-capacity
C_epsilon = C_s(loop*error)/log(2);

for n=nn
    error_rate = 0;
    rate = rate_previous-99/loop1;
    while (error_rate < error)
        %Normal Approximation for the error at a given rate
        error_rate = mean(1 - qfunc( sqrt(n)*(rate-C_array)./sqrt(V_array) ));
        rate = rate + 1/loop1;
    end
    
    rate_previous=rate;

    rate_na=[rate_na,rate/log(2)];
end