function rate_na_lt = awgn_na_lt(snr_db,nn,error)


P=10^(snr_db/10);

rate_na_lt=zeros(1,length(nn));

%find the minimum blocklength required for the long-term power constraint to be beneficial  
C=log(1+P);
V = 1-1/(1+P)^2;
n_min = ( (1+P)/P*sqrt(2*3.1415926*V)*(1-error)*exp(qfuncinv(error)^2/2) + qfuncinv(error)/(1+P)^2/sqrt(V))^2;

n_index=find(nn <= n_min);

% for every n<=n_min, the normal approx. for the st power constraint coincides 
if length(n_index)>0
    rat_na_lt(n_index) = (log(1+P) - sqrt(1-(1+P)^(-2))./sqrt(nn(n_index)) + 0.5*log(nn(n_index))./nn(n_index))/log(2); 
end
    


R_st = C - sqrt(V./nn)*qfuncinv(error);
C_lt=log(1+P/(1-error)); 

index_large_n = find(nn>n_min);

for ii=index_large_n
    n=nn(ii);
    for R_star = R_st(ii):0.0001:C_lt
        omega = P: 0.00001:P*1.2; %note that larger endpoing for omega (P*1.2) gives a better estimate for rate_na;
        
        q_arg=sqrt(n)*(log(1+omega)-R_star)./sqrt(1-(1+omega).^(-2)) ;
        
        error_temp = min( 1-P./omega + P.*qfunc( q_arg)./omega);
        if error_temp> error
            rate_na_lt(ii)= (R_star +0.5*log(n)./n)./log(2) ;
            break
        end
    end
end
