%rate_c = converse_simo(nn,P,error,rx,K) computes the converse bound on the
%maximal channel coding rate over a quasi-static single-input
%multiple-output (SIMO) Rician fading channel with full channel state information.  
%
%The function takes the following inputs:
%
%nn: blocklength, scalar/vector;
%P: input power (in linear scale), scalar;
%error: target block error probability, scalar, should be greater than 10^(-6);
%rx: number of receive antennas, scalar;
%K: rician K-factor, scalar; the default value is 0
function rate_c = converse_simo(nn,P,error,rx,K)

if (nargin < 5) || isempty(K)
	K = 0;
end

loop = min(1000/error, 10^7);    % number of samples of channel gain for Monte Carlo

%Generate samples of channel gain for Monte Carlo simulation
%channel gain is G = |H_1|^2 + ... +|H_{rx}|^2
if rx ==1
    G= abs(sqrt(K/(K+1))+sqrt(1/(K+1)/2)*(randn(rx,loop)+1i*randn(rx,loop))).^2 ;
else
    G = sum(abs(sqrt(K/(K+1))+sqrt(1/(K+1)/2)*(randn(rx,loop)+1i*randn(rx,loop))).^2 );
end

rate_c=[];

loop21=1000; %number of samples in numerical integration, set larger value for better precision
loop22=1000;


%samples of channel gain for numerical computation
%Todo: optimize the way of sampling G from ncx2 (noncentral chi square)
g_mid = ncx2inv(error*2, 2*rx,2*rx*K)/(2*K+2);
g_max = ncx2inv(1-10^(-5),2*rx,2*rx*K )/(2*K+2);

G1 = (g_mid/loop21):(g_mid/loop21):g_mid;
G2 = (g_mid+ (g_max -g_mid)/loop22):(g_max -g_mid)/loop22:g_max;
G_num=[G1,G2];

%compute the pdf of channel gain
pdf_G = ncx2pdf(G_num*(2*K+2), 2*rx, 2*rx*K)*(2*K+2);

%pdf_G = (K+1) * exp(sqrt(4*K*2*(K+1)*G)-(K+1).*G-2*K).* sqrt(2*(K+1).*G/4/K).*besseli(1,sqrt(4*K*2*(K+1)*G),1);

precision = 10000; %number of samples for numerical integration; set larger value to increase precision

for n = nn
    
    %compute gamma_n through Monte Carlo
    S_n = n.*log(1+P*G) + n - ncx2rnd(2*n, 2*n/P./G).*P.*G./(1+P*G)/2; %information density
    S_n_sort=sort(S_n);
    gamma_n = S_n_sort(error*loop); % Pr[S_n <= gamma_n]=error
    %%%Todo: Change to numerical computations
    
    % Calculate Q[log dP/dQ >=gamma_n]; Given G, log dP/dQ is ncx2 distributed
    y_a=zeros(1,length(G_num)); % conditional probability P(Ln\geq n\gamma_n) given g for different g values
     
    for ii=1:length(G_num)
        gg= G_num(ii);
        k_n = 2/P/gg*( n*log(1+P*gg) + n-gamma_n );
        noncentrality = 2* n*(1+P*gg)/P/gg;
        %use ncx2cdf to compute cdf of noncentral chi-square
        cond_cdf = ncx2cdf(k_n , 2*n, noncentrality);
        
        if cond_cdf ==0 && k_n>0
            %Numerically compute the cdf of ncx2; this is done by convolving chi^2 pdf and |Gaussian(noncentrality,1)|^2
            
            step = sqrt(k_n)/precision;
            
            t = (sqrt(noncentrality)-sqrt(k_n)) :step:(sqrt(noncentrality) + sqrt(k_n));
            
            y_a(ii) = sum( exp(-1*(t).^2/2 + log( gammainc((k_n -(sqrt(noncentrality) - t).^2 )/2,(2*n-1)/2) ) ) ) * step./sqrt(2*pi); 
        else
            y_a(ii)=cond_cdf;
        end
    end
    
    beta = sum( y_a(1:length(G1)) .* pdf_G(1:length(G1)))*g_mid/loop21 + sum(y_a(length(G1)+1: end) .*pdf_G(length(G1)+1:end)) *(g_max- g_mid)/loop22; % average over G.
    
    rate=-1*log2(beta)/(n-1);
    rate_c=[rate_c,rate];
end