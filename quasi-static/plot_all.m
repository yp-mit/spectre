function [C_e, rate_a,rate_c,rate_na] = plot_all(nn,P,error,tx,rx,K)
%[C_epsilon, rate_a,rate_c,rate_na] = plot_all(nn,P,error,rx,K) plots
%converse/achievability/normal approximation/capacity plots on the same
%figure.  
%
%The function takes the following inputs:
%
%nn: blocklength, scalar/vector;
%P: input power (in linear scale), scalar;
%error: target block error probability, scalar, should be greater than 10^(-6);
%tx: number of transmit antennas, scalar;
%rx: number of receive antennas, scalar;
%K: rician K-factor, scalar; the default value is 0

if (nargin < 6) || isempty(K)
	K = 0;
end

if tx==1
    cd ./SIMO_rician;
    
    g0 = ncx2inv(error,2*rx, 2*rx*K)/(2*K+2); %compute the error-th quantile of the sum of rx rician distributions with K-factor K
%
    C_e = log2(1 + g0*P);% epsilon-capacity
    
    rate_c = converse_simo(nn,P,error,rx,K);
    rate_a1= ach_simo_nocsi(nn,P,error,rx,K);
    rate_a2= ach_simo_csir(nn,P,error,rx,K);
    rate_na=normapprox_simo(nn,P,error,rx,K);
    rate_a = [rate_a1;rate_a2];
    
    plot(nn, C_e*ones(size(nn)), 'k-*');
    
    hold on;
    grid on;
    plot(nn, rate_c,'r');
    plot(nn,rate_a1,'b--');
    plot(nn,rate_a2,'b');
    plot(nn,rate_na,'k');
    
    legend(['Capacity = ',num2str(C_e)],'Converse','Achievability (no CSI)', 'Achievability (CSIR)' ,'Normal approximation');
else
    cd ./MIMO_ISO;
    rate_a = mimo_iso_ach(nn,P,error,tx,rx,K);
    rate_c = mimo_iso_conv(nn,P,error,tx,rx,K);
    [rate_na,C_e]= normapprox_mimo_iso(nn,P,error,tx,rx,K);
    
    plot(nn, C_e * ones(size(nn)),'k-*');
    hold on 
    grid on
    plot(nn,rate_c,'r');
    plot(nn,rate_a,'b');
    plot(nn,rate_na,'k');
    legend(['Capacity = ', num2str(C_e)],'Converse','Achievability (no CSI)','Normal approximation');
end

xlabel('Blocklen, n'); ylabel('Rate, bits/(ch. use)');
cd ..;