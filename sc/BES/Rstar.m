function out = Rstar(d, n, e, alpha, fun)
%finds excess distortion for a given rate at a given blocklength for erased fair coin
%flip source
%d - excess distortion
%n - block length (scalar)
%e - excess probability
%alpha - erasure rate
%fun - which function to use for calculation

%
%   Created in 2012 by Victoria Kostina (vkostina@caltech.edu)
%

%starting points - 'persistent' to make optimization faster
persistent x0;

Rd = RateDistortion(alpha, d);
accuracyC = 1;
slack = 0.5;

options = optimset('Display', 'off');
switch lower(fun)
    case 'spherecoveringa'
        if isempty(x0) 
            x0 = [Rd 5]; %rate
        end
        out = SphereCoveringA();
    case 'spherecoveringc'
        %converse via sphere covering
        if isempty(x0) 
            x0 = [Rd 5]; %rate
        end
        out = SphereCoveringC();
    case 'normal'
        %normal approximation
        out = Normal();
    case 'aeq'
        if isempty(x0)
            x0 = [Rd 5]; %rate
        end
        out = Aeq();
    case 'ceq'
        %converse via sphere covering
        if isempty(x0)
            x0 = [Rd 5]; %rate
        end
        out = Ceq();
    case 'normaleq'
        %normal approximation
        out = Normaleq();
        
    otherwise
        disp('Error: unknown type.')
end



%--------------------------------------------------------------------------
    function out = SphereCoveringA()
        %achievability
        
        out = fzero(@(x)Peub(x,n,d) - e, x0, options);
        updatex0(out);
        function out = Peub(R, n, d)
            %excess probability
            M = 2^(n*R);
            out = 0;
            for t = 0:n
                out = out + binopdf(t,n,alpha)*fail(n,t);
            end
            %out = fail(n,alpha*n);
            fprintf('SphereCoveringA: n = %i, R = %f, Peub = %f\n', n, R, out);
            
            function out = fail(n,t)
                %Probability that distortion exceeds d
                out = 0;
                for a = 0:t
                    out = out + C(t, a, 'ub')*nonerased(floor(d*n) - a);
                end
                out = out/2^t;
                
                function out = nonerased(b)
                    % Probability that distortion in nonerased symbols exceeds
                    % d - a
                    out = 2^(-n + t)*Csum(n - t, b, 'lb', accuracyC);
                    if n <= 50
                        out = (1 - out)^M;
                    else
                        out = exp(-M*out);
                    end
                    
                end
            end
        end
    end

%--------------------------------------------------------------------------
    function out = SphereCoveringC()
        %converse
        
        out = fzero(@(x)Pelb(x,n,d) - e, x0, options);
        updatex0(out);
        function out = Pelb(R, n, d)
            %excess probability
            out = 0;
            for t = 0:n
                out = out + binopdf(t,n,alpha)*success(n,t);
            end
            out = 1 - out;
            fprintf('SphereCoveringC: n = %i, R = %f, Peub = %f\n', n, R, out);
            
            function out = success(n,t)
                %number of sequences with distortion within d
                out = 0;
                for a = 0: floor(min(d*n, t))
                    nonerased = Csum(n - t,floor(d*n) - a, 'ub', accuracyC)*2^(n*R - n + t);
                    if nonerased > 1
                        nonerased = 1;
                    end
                    out = out + C(t, a, 'ub')*nonerased;
                end
                out = out/2^t;
                %fprintf('SphereCoveringC-success: n = %i, e = %f, P = %f\n', n, t/n, out*2^(n*(R - 1)));
            end
        end
    end

%--------------------------------------------------------------------------
    function out = Normal()
        %normal approximation
        lambda = log2( (d - alpha/2)/(1 - alpha/2 - d));
        V = alpha*(1 - alpha)*(log2( (2^(-lambda/2) + 2^(lambda/2))/2 ) )^2 + alpha/4*lambda^2;
        out = RateDistortion(alpha, d) + sqrt(V/n)*Qinv(e);
    end


%--------------------------------------------------------------------------
    function out = Aeq()
        %achievability for the asymptotically equivalent problem
        
        out = fzero(@(x)Peub(x,n,d) - e, x0, options);
        updatex0(out);
        function out = Peub(R, n, d)
            %excess probability
            M = 2^(n*R);
            out = 1 - binocdf(floor(2*n*d), n, alpha);
            for t = 0:floor(2*n*d)
                out = out + binopdf(t,n,alpha)*fail(n,t);
            end
            %out = fail(n,alpha*n);
            %fprintf('SphereCoveringA: n = %i, R = %f, Peub = %f\n', n, R, out);
            
            function out = fail(n,t)
                %Probability that distortion exceeds d
                out =  nonerased(floor(d*n - t/2) );
                
                function out = nonerased(b)
                    % Probability that distortion in nonerased symbols exceeds
                    % d - a
                    out = 2^(-n + t)*Csum(n - t, b, 'lb', accuracyC);
                    if n <= 50
                        out = (1 - out)^M;
                    else
                        out = exp(-M*out);
                    end
                    
                end
            end
        end
    end

%--------------------------------------------------------------------------
    function out = Ceq()
        %converse for the asymptotically equivalent problem
        
        out = fzero(@(x)Pelb(x,n,d) - e, x0, options);
        updatex0(out);
        function out = Pelb(R, n, d)
            %excess probability
            out = 0;
            for t = 0:floor(2*n*d)
                out = out + binopdf(t,n,alpha)*success(n,t);
            end
            out = 1 - out;
            %fprintf('SphereCoveringC: n = %i, R = %f, Peub = %f\n', n, R, out);
            
            function out = success(n,t)
                %number of sequences with distortion within d
                out = min(1,  Csum(n - t,floor(d*n - t/2), 'ub', accuracyC)*2^(n*R - n + t));
                %fprintf('SphereCoveringC-success: n = %i, e = %f, P = %f\n', n, t/n, out*2^(n*(R - 1)));
            end
        end
    end



%--------------------------------------------------------------------------
    function out = Normaleq()
        %normal approximation
        lambda = log2( (d - alpha/2)/(1 - alpha/2 - d));
        V = alpha*(1 - alpha)*(log2( (2^(-lambda/2) + 2^(lambda/2))/2 ) )^2;
        out = RateDistortion(alpha, d) + sqrt(V/n)*Qinv(e);
    end


%--------------------------------------------------------------------------
%Utilities for faster computation
    function updatex0(x)
        if ~isnan(x)
            x0(2) = x + slack;
        end
    end

end
