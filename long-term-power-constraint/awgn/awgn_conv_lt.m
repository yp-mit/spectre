%converse bound for AWGN channel under a long-term power constraint
function [rate_conv_lt]=awgn_conv_lt(snr_db, nn, error)


rate_na_lt =  awgn_na_lt(snr_db,nn,error)*log(2); %normal approximation 

Hb_error = -error*log(error) -(1-error)*log(1-error);
P=10^(snr_db/10);


C=log(1+P);
V=1-1/(1+P)^2;

rate_conv_lt=zeros(1,length(nn));
 
for index_n=1:length(nn)
    
    n = nn(index_n);
    disp(['awgn_conv_lt(): n=',num2str(n)]);
     
    error_array=[];
    

    gamma_array=rate_na_lt(index_n):0.0001:C*1.1;
    
    for gamma=gamma_array
        gamma
        xx = exp(gamma)-1:0.0001:exp(gamma)-0.8; %estimate the range of the tangent pt. Modify this range according to the parameters n, P, and error;

       y=marcumq(sqrt(2*n./xx), sqrt(n*(log(1+xx)+1 - gamma).*(1+xx)*2./xx  ), n);
       
       slope = (1-y)./xx;
       index = find (slope == max(slope));
       x0 = xx(index(1)); %compute the tangent pt
       
       error_n = 1-P/x0 + P/x0*y(index(1));
       
       error_array=[error_array,error_n];
       
       if error_n> 1.5*error %this number 1.5 should be changed according to the parameters n, P, and error to avoid unnecessary computations. 
           break;
       end
    end

    index_error= find(error_array > error);   
    varrho_n = log(1+sqrt(n/2/pi)*( log(P) + (n*C + Hb_error)/(1-error) + log(1+exp(-(n*C + Hb_error)/(1-error) )/P)) )/n;
    
    rate_array = gamma_array(index_error) - log(error_array(index_error)-error)/n + varrho_n;
    
    rate_conv_lt(index_n) = min(rate_array)./log(2);     
end
