%rate_a=ach_simo_csir(nn,P,error,rx,K) computes the achievability bound on
%the maximal channel coding rate over a single-input multiple-output (SIMO)
%quasi-static Rician fading channel with perfect channel state information at
%the receiver (CSIR).
%
%The function takes the following inputs:
%
%nn: blocklength, scalar/vector;
%P: input power (in linear scale), scalar;
%error: target block error probability, scalar, should be greater than 10^(-6);
%rx: number of receive antennas, scalar;
%K: rician K-factor, scalar; the default value is 0.
function rate_a = ach_simo_csir(nn,P,error,rx,K)

if (nargin < 5) || isempty(K)
	K = 0;
end

loop =  min(1000/error, 10^7);
loop =  max(loop, 10^5);

if rx ==1
    G= abs(sqrt(K/(K+1))+sqrt(1/(K+1)/2)*(randn(rx,loop)+1i*randn(rx,loop))).^2 ;
else
    G = sum(abs(sqrt(K/(K+1))+sqrt(1/(K+1)/2)*(randn(rx,loop)+1i*randn(rx,loop))).^2 );
end


rate_a=[];

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


for n=nn
    
    %compute gamma_n using Monte Carlo;
    %Todo: use numerical methods
    S_n = n.*log(1+P*G) + n - ncx2rnd(2*n, 2*n/P./G).*P.*G./(1+P*G)/2;
    S_n_sort=sort(S_n);
    tau=2*error/n;
    gamma_n = S_n_sort(floor((error-tau)*loop));
    %%%Todo1: Change to numerical computations
    %%%Todo2: optimize over tau
    
    
    %Next we lower-bound kappa_tau
    %According to [MolavianJazi&Laneman13], dP/dQ is upper bounded by a universal constant ([MolavianJazi&Laneman13] deals with AWGN, but it is shown there that the bound is independent of the input power)
    y = ncx2rnd(2*n, 2 * n * P*G)/2; %noncentral chisquare dist. with 2n dof and noncentral para. 2npG
    log_dpdq = n*log(1+P*G) -(sqrt(y) -sqrt(n*P*G)).^2 + y./(1+P*G) - (n-1)/2*log(n.*G.*P.*y) + gammaln(n) + log(besseli(n-1, 2*sqrt(n.*G*P.*y),1));
    
    %%use the following upper bound for log_dpdq whenever n is large, and besseli function gets overflow
    %dp_a = log(sqrt(pi)/4) + (1/4 - n/2) * log( n*P.*y.*G) + n *log(1+P.*G) + gammaln(n) - 0.25*log(1+(n-1)^2./(4*n*P.*G.*y)) - n.*G*P - P.*G./(1+P.*G).* y - asinh((n-1)./sqrt(4*n*P.*G.*y))*(n-1) + 2*sqrt(n*P*G.*y).*sqrt(1+(n-1)^2./(4*n*P*G.*y));
    
    max_ldpdq=max(log_dpdq); %maximal of log (dP/dQ)
    %Todo: actually max_ldpdq can be upper-bounded analytically for each n. 
    
    
    
    precision=10000;
    %next, we compute beta
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
            
            y_a(ii) = sum( exp(-1*(t).^2/2 + log(gammainc((k_n -(sqrt(noncentrality) - t).^2 )/2,(2*n-1)/2)))) * step./sqrt(2*pi);
        else
            y_a(ii)=cond_cdf;
        end
    end
    beta = sum( y_a(1:length(G1)) .* pdf_G(1:length(G1)))*g_mid/loop21 + sum(y_a(length(G1)+1: end) .*pdf_G(length(G1)+1:end)) *(g_max- g_mid)/loop22; % average over G.
    
    %kappa_tau >= 1/gamma (tau - P[log dP/dQ >=log gamma]) where P[log dP/dQ >=log gamma]=0 if log gamma = max_ldpdq.
    rate=(log2(tau)-log2(beta) -max_ldpdq./log(2))/n;
    rate_a=[rate_a,rate];
    
end
