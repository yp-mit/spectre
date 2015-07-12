%info_bits=function conv_epb_csir(EE, epsil);
%
%This function computes the converse bound on the maixmal number of
%information bits that can be transmitted over a Rayleigh-fading channel
%(perfect CSIR) with a fixed energy budget "EE" and error probability
%"epsil". Here, "EE" is a vector, and "epsil" is a scalar.

function info_bits=conv_epb_csir(EE, epsil)


info_bits=[];
for ij = 1:length(EE)
    
    E = EE(ij);

    eta_array =[];
    pr_array=[];
    
    for eta =10.^(log10(E)-2:0.01: log10(E)+1)
        
        if eta<=6
            error_t =  qfunc((E-eta)./sqrt(2*E));
            
            if error_t >epsil
                eta_array =[eta_array, eta];
                pr_array=[pr_array, error_t];
            end
            
        else
            
            x_array=eta:0.01:10*eta;
            y_array = qfunc((eta -x_array)./sqrt(2*x_array)) - 1/4/sqrt(pi)*exp(-(x_array-eta).^2./x_array/4) .*(eta./sqrt(x_array) + sqrt(x_array));
            index_tangent = length( find(y_array < 0 ) );
            if index_tangent == 0
                break;
            end
            
            x_0 = x_array(index_tangent);
            
            if x_0 < E
                error_t = qfunc((E-eta)./sqrt(2*E));
            else
                error_t= 1- qfunc((eta-x_0)/sqrt(2*x_0))*E/x_0;
            end
            
            if error_t >epsil
                eta_array =[eta_array, eta];
                pr_array=[pr_array, error_t];
            end
            
        end
    end
    
    logM = min(eta_array - log(pr_array -epsil));
    info_bits =[info_bits, logM/log(2)];
    
%     Eb=log10(E/logM*log(2))*10;
%     disp(['k^*(E,error)=',num2str(logM/log(2)),';    Eb=',num2str(Eb)]);
end

