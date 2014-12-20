%rate_a = mimo_iso_ach(nn,P,error,tx,rx,K) evaluates the achievability
%bound on the maximal achievable rate over a quasi-static MIMO Rician
%fading channel with no CSI and with isotropic inputs.  
%
%The function takes the following inputs:
%
%nn: blocklength, scalar/vector;
%P: total input power (linear scale), scalar;
%error: block error probability, scalar;
%tx: number of transmit antennas, scalar;
%rx: number of receive antennas, scalar;
%K:  rician K factor, scalar; the default value is 0.
function rate_a = mimo_iso_ach(nn,P,error,tx,rx,K)

if (nargin < 6) || isempty(K)
	K = 0;
end

loop_h = round(1000/error); %number of Monte Carlo samples

P=P/tx; %transmit power per transmit antenna

tau = [10:(error*loop_h-10)]./loop_h;

rate_a=[];

for  n= nn
    
    sin_arr=zeros(1,loop_h);
    
    for ii=1:loop_h
        h=sqrt(0.5/(K+1))*(randn(tx,rx)+1i*randn(tx,rx)) + sqrt(K/(K+1));
        %X_0 = [sqrt(n*P)*eye(tx),0];
        y1 =sqrt(n*P)*eye(tx)* h  +  sqrt(0.5)*(randn(tx,rx) + 1i*randn(tx,rx));
        %
        W2=sqrt(0.5)*(randn( n-tx,rx) + 1i*randn(n-tx,rx));
        QAQB=(y1'*y1) * inv(y1'*y1+W2'*W2);
        sin_arr(ii)=abs(det(eye(rx) - QAQB));
    end
    
    sin_sort=sort(sin_arr,'descend');
    
    logM_n=[];
    %chernoff bound (optimizing over alpha)
    for alpha=2 : n-rx-tx
        logM = max( log(tau) - alpha*log(sin_sort(floor((error-tau)*loop_h))) );
        for index_rx =1:rx
            logM = logM - sum(log( n-index_rx-[0:1:tx-1] )) + sum(log(n -alpha -index_rx - [0:1:tx-1]));
            logM_n=[logM_n,logM];
        end
        
    end
    rate = max(logM_n)/n/log(2);
    rate_a =[rate_a,rate];

end
