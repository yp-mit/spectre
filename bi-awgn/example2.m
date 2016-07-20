%% 2) PACKET ERROR RATE BOUNDS
disp(' ')
disp('PACKET ERROR RATE BOUNDS')

n = [5e2 1e3 2e3 5e3 1e4 1e5 1e6];
Pe = 10.^([-0.25:-0.25:-1.75,-2:-0.5:-7]);
uno = ones(size(Pe));
R = 1/2;
tic

if (0) % set to 1 if you want to generate values
    
    % NA bound 
    % expected execution time 110 sec (2 min)
    for k = 1:length(n)    
        Om_NA(k,:) = biawgnPPVbound(n(k)*uno,Pe,R*uno,'normal','error');
    end; toc; tic
    
    % PPV bound 
    % expected execution time 721 sec (12 min)
    for k = 1:length(n)    
        Om_PPV(k,:) = biawgnPPVbound(n(k)*uno,Pe,R*uno,'On2','error');
    end; toc; tic

    save('example2.mat')
else
    load('example2.mat')
end

figure(1)
subplot(1,1,1)
set(0,'defaulttextinterpreter','latex')
semilogy(Om_NA,Pe,'--')
hold on
semilogy(Om_PPV,Pe)
hold off
xlabel('SNR $\Gamma$')
ylabel('error probability $P_e$')
grid