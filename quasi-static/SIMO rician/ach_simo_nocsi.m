%rate_a=ach_simo_nocsi(nn,P,error,rx,K) computes the achievability bound on
%the maximal channel coding rate over a single-input multiple-output (SIMO)
%quasi-static Rician fading channel with no channel state information at
%transmitter or receiver.
%
%The function takes the following inputs:
%
%nn: blocklength, scalar/vector;
%P: input power (in linear scale), scalar;
%error: target block error probability, scalar, should be greater than 10^(-6);
%rx: number of receive antennas, scalar;
%K: rician K-factor, scalar; the default value is 0.
function rate_a=ach_simo_nocsi(nn,P,error,rx,K)

if (nargin < 5) || isempty(K)
	K = 0;
end

loop=min(1000/error, 10^7);

rate_a=[];
hh = sqrt(0.5/(K+1))*(randn(rx,loop)+1i*randn(rx,loop)) + sqrt(K/(K+1)); %samples of Rician fading channels
w1 = sqrt(0.5) * (randn(rx, loop) + 1i*randn(rx,loop)); %AWGN noise

for n=nn
    data_a=zeros(1,loop);
    
    %By spherical symmetry, it suffices to consider x0=[(nP)^{-1/2},0,...,0]
    y1 = hh.*sqrt(n*P) + w1; %The received vector during the first channel use
    for ii=1:loop
        y1_sample = y1(:,ii);
        W = sqrt(0.5) * (randn(rx, n-1) + 1i*randn(rx,n-1));

        %projection of x0 onto Y
        data_a(ii) = abs(y1_sample'* inv(y1_sample*y1_sample' + W*W')*y1_sample);
        %
    end

    
    data_sort = sort(data_a); 
    
    
    %maximize over tau...
    tau_array = [floor(error*loop)-1:-1:1]/loop ;
    
    gamma_tau = 1 - data_sort(1:floor(error*loop)-1);
    
    log_Fn = log(betainc(gamma_tau, n-rx,rx));  
    index_set= find (log_Fn < -10^100);%detect if the output of betainc is zero
    
    if length(index_set)>0
       log_Fn(index_set) = (n-rx) * log(gamma_tau(index_set))  - sum(log(n-rx+1:1:n-1)) + log(gamma(rx)) ;
    end
    
    %
    log_M = max(log(tau_array) - log_Fn );
    
    rate_a =[rate_a, log_M/n/log(2)];
end
