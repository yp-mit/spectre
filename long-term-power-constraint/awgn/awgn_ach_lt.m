%achievability bound for AWGN channel under a long-term power constraint
%

function rate_ach_lt=awgn_ach_lt(snr_db, nn, error)

snr = 10^(snr_db);

rate_ach_lt = zeros(1,length(nn)); %achievable rate (kappa beta)

%Check whether a long-term power constraint will be beneficial 
C=log(1+snr);
V = 1-1/(1+snr)^2;
n_min = ( (1+snr)/snr*sqrt(2*pi*snr)*(1-error)*exp(qfuncinv(error)^2/2) + qfuncinv(error)/(1+snr)^2/sqrt(V))^2;

n_index=find(nn <= n_min);

if length(n_index)>0
    rate_ach_lt(n_index) = kappabeta_ach(2*nn(n_index),error,snr,1)./nn(n_index); %use the kappa beta bound for short-term power constraint  
end


index_large_n = find(nn>n_min);

number_taus=20; %increase this number of gain tightness (at the cost of complexity)
number_errors=20; %increase this number of gain tightness (at the cost of complexity)

for index_n= index_large_n
    n=nn(index_n);
    %tic
    disp(['ach_awgn_lt():  n=', num2str(n)]);

    temp_rate_arr=[];
    for error_n = (error/number_errors): (error /number_errors): (error * (number_errors-1)/number_errors)
        %ach SISO kappa beta bound
        P_n= snr*(1-error_n)./(1-error);
        tau_array = (error_n/number_taus) : (error_n/number_taus): (error_n*(number_taus-1)/number_taus);
        for tau = tau_array
            temp_rate = (log2(kappa_inf(tau, P_n)) - betaq_up_v2(1-error_n+tau, 2*n, P_n))/n;
            temp_rate_arr=[temp_rate_arr,temp_rate];
        end
        rate_ach_lt(index_n)=max(temp_rate_arr); 
    end
end