function out = Rstar(p, e, n, fun)
%finds average rate compatible with given error probability at a given blocklength for coin
%flip source
%p - source bias
%e - excess probability
%n - block length (scalar)
%fun - which function to use for calculation

%
%   Created in 2013 by Victoria Kostina (vkostina@caltech.edu)
%



R0 = h(p);
V = p*(1-p)*(log2((1-p)/p))^2;



switch lower(fun)
    case 'approx'
        %gaussian approximation
        out = Approx();
    case 'exact'
        %Exact average rate 
        out = Exact();
    case 'exacti'
        %exact min I(S; Z) subject to P[S neq Z] \leq e
        out = ExactI(); 
    case 'ie'
        %E[i_S(S, epsilon)]
        out = ie();
    otherwise
        disp('Error: unknown type.')
end        







%--------------------------------------------------------------------------
    function out = Approx()
        %Gaussian approximation
        out = (1 - e)*R0 - sqrt(V/n/2/pi)*exp( - (Qinv(e))^2/2); % - (1-e)/2*log2(n)/n;
    end





%--------------------------------------------------------------------------
    function out = Exact()
       %Exact average rate 
        
        M = 0; %the number of previous sequences
        j = 0; %j is the target length
        Nj = 0;%total number of sequences of length j
        Elen = 0;
        cdf = 0; %total probability
        
        for k = 0:n
            m = C(n, k, 'ub'); %m is the number of remaining sequences of current type
            prob = 2^(k*log2(p) + (n - k)*log2(1-p));
            
            while m > 0
                if floor(log2(M + m)) > j
                    %2^j - Nj sequences of type k will have length j
                    dm = 2^j - Nj;
                else
                    %all remaining sequences of type k will have length j
                    dm = m;
                end
                %dm is the number of sequences of type k and length j
                %fprintf('%i seq of type %i with length %i \n', dm, k, j );
                
                pdf = 2^(log2(dm) + log2(prob));
                cdf = cdf + pdf;
                
                if cdf > 1 - e
                    w = (1 - e - cdf + pdf)/pdf; %encoder randomization
                    Elen = Elen + j*w*pdf;
                    out = Elen / n;
                    return;
                end
                
                Elen = Elen + j*pdf;
                
                m = m - dm;
                M = M + dm;
                Nj = Nj + dm;
                
                
                if log2(Nj) == j
                    % exhausted all codewords of length j
                    j = j + 1;
                    Nj = 0;
                end
            end
            
        end
        
        %mysum = Csum(n, rstar, 'ub', 1)
        
        out = Elen / n;
        
    end










%--------------------------------------------------------------------------
    function out = ExactI()
        %exact min I(S; Z) subject to P[S neq Z] \leq e
        
        M = 0; %the number of masses with probability strictly greater than prob
        Ei = 0;
        cdf = 0; %total probability
        eta = 1;
        etastar = e/(2^n-1);
        % 
        
        for k = 0:n
            m = C(n, k, 'ub'); %m is the number of sequences of current type
            prob = 2^(k*log2(p) + (n - k)*log2(1-p));
   
            
            if cdf - 2^(log2(max(0, M-1)) + log2(eta)) > 1 - e
               etastar = 2^( log2(cdf - 1 + e) - log2(M-1) ); % etastar < eta
               break;
            end

            eta = prob;
            M = M + m;
            
            pdf = 2^(log2(m) + log2(prob));
            cdf = cdf + pdf; 
            
            Ei = Ei - pdf*log2(prob);
                      
        end
        
        out = ( Ei - (1 - cdf)*log2(etastar) + e*log2(etastar)+ (1 - e)*log2(1 - e) ) / n;
        
    end











%--------------------------------------------------------------------------
    function out = ie()
        %E[i_S(S, epsilon)]        

        Elen = 0;
        cdf = 0;
        
        for k = 0:n
            m = C(n, k, 'ub'); %m is the number of sequences of current type
            prob = 2^(k*log2(p) + (n - k)*log2(1-p));
            pdf = 2^(log2(m) + log2(prob));
            cdf = cdf + pdf;    
                if cdf > 1 - e
                    w = (1 - e - cdf + pdf)/pdf; %encoder randomization
                    Elen = Elen - log2(prob)*w*pdf;
                    out = Elen / n;
                    return;
                end
                
                Elen = Elen - log2(prob)*pdf;                         
        end    
        
        out = Elen / n;
        
    end


end




