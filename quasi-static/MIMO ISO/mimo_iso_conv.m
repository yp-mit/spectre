%rate = mimo_iso_conv(nn,P,error,tx,rx,K) computes the converse bound on
%the maximal channel coding rate over a quasi-static MIMO Rician fading
%channel with isotropic inputs. 
%
%The function takes the following inputs:
%
%nn: blocklength, scalar/vector;
%P: total input power (linear scale), scalar;
%error: block error probability, scalar;
%tx: number of transmit antennas, scalar;
%rx: number of receive antennas, scalar;
%K:  rician K factor, scalar; the default value is 0.

function rate_c = mimo_iso_conv(nn,P,error,tx,rx,K)
if (nargin < 6) || isempty(K)
	K = 0;
end
loop_h=round(1000/error);

P=P/tx;

rate_c=[];

min_tx_rx = min(tx,rx);

eig_vec=zeros( min_tx_rx, loop_h);

for ij = 1:loop_h
    h = (randn(tx,rx)+1i*randn(tx,rx))*sqrt(0.5/(K+1)) + sqrt(K/(K+1));
    %eigenvalues of H^{H}*Q*H, where Q is a scaled identity matrix P*I
    eig_vec(:,ij)= P* abs(svd( h)).^2;
end

%nn=[60:10:200 220:20:500 550:50:1000];

for n=nn
    %information density under P
    S_n = sum( n*log(1 + eig_vec)) + n*min_tx_rx  - sum(ncx2rnd(2*n, 2*n ./ eig_vec).* eig_vec./ (1 + eig_vec))/2;
    
    S_n_sort=sort(S_n);

    tau = (error + 1/loop_h ): 1/loop_h : 2*error;
    
    gamma_n_array = S_n_sort(round(tau*loop_h));
    
    %lower-bounding beta_{error} using the standard inequality 
    y_a = gamma_n_array/n - log(tau-error)/n;
    %optimizing over gamma
    y_min=min(y_a);
    
    rate =y_min./log(2);
    rate_c=[rate_c,rate];
end



%