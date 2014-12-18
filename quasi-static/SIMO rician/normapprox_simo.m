%rate_na = normapprox_simo(nn,P,error,rx,K) computes the normal
%approximation of the maximal channel coding rate over a quasi-static SIMO Rician-fading channel. 
%
%The function takes the following inputs:
%
%nn: blocklength, scalar/vector, should be put in an ascending order;
%P: input power (in linear scale), scalar;
%error: target block error probability, scalar;
%rx: number of receive antennas, scalar;
%K: rician K-factor, scalar; the default value is 0.
function rate_na = normapprox_simo(nn,P,error,rx,K)

if (nargin < 5) || isempty(K)
	K = 0;
end

g0 = ncx2inv(error,2*rx, 2*rx*K)/(2*K+2); %compute the error-th quantile of the sum of rx rician distributions with K-factor K
%
C_error = log(1 + g0*P);% epsilon-capacity

rate_na=[];

loop1=10000;

loop21=5000; %number of samples in numerical integration, set larger value for better precision
loop22=5000;

%samples of channel gain for numerical computation
%Todo: optimize the way of sampling G from ncx2 (noncentral chi square)
g_mid = ncx2inv(error*2, 2*rx,2*rx*K)/(2*K+2);
g_max = ncx2inv(1-10^(-5),2*rx,2*rx*K )/(2*K+2);

G1 = (g_mid/loop21):(g_mid/loop21):g_mid;
G2 = (g_mid+ (g_max -g_mid)/loop22):(g_max -g_mid)/loop22:g_max;
G_num=[G1,G2];

%compute the pdf of channel gain
pdf_G = ncx2pdf(G_num*(2*K+2), 2*rx, 2*rx*K)*(2*K+2);

error=0.001;

rate_previous=100/loop1;
for n=nn
    for rate=rate_previous-99/loop1 : 1/loop1 :100
        phi_G = 1 - qfunc(sqrt(n)*(rate-log(1+G_num*P))./sqrt(1-1./(1+P*G_num).^2));
        error_rate = sum(phi_G(1:length(G1)) .* pdf_G(1:length(G1)))*g_mid/loop21 + sum(phi_G(length(G1)+1: end) .*pdf_G(length(G1)+1:end)) *(g_max- g_mid)/loop22; % average over G.
        if (error_rate>=error)
            break;
        end
    end
    rate_previous=rate;
    rate2 = (rate+1/2/n*log(n))./log(2);
    rate_na=[rate_na,rate2];
    %save rate_a rate_a nn
end