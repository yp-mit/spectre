if(0) % set to (1) if you wish to run the biawgnPPVbound function
    % SPECTRAL EFFICIENCY CURVES - UPPER BOUNDS ON RATE

    n = 1e3;
    Pe = 1e-5;
    SNRdB = -15:10;
    uno = ones(size(SNRdB));tic

    % normal approximation - expected execution time 2 sec
    rhoNA = 2*biawgnPPVbound(n*uno,Pe*uno,SNRdB,'normal'); % spectral efficiency 
    EbN0NA = SNRdB-10*log10(rhoNA);toc; tic
    % O(n^-2) approximation - expected execution time 90 sec
    rhoPPV = 2*biawgnPPVbound(n*uno,Pe*uno,SNRdB,'On2'); % spectral efficiency 
    EbN0PPV = SNRdB-10*log10(rhoPPV);toc; tic
    % O(n^-3) approximation - expected execution time 120 sec
    rhoPPVb = 2*biawgnPPVbound(n*uno,Pe*uno,SNRdB,'On3'); % spectral efficiency 
    EbN0PPVb = SNRdB-10*log10(rhoPPVb);toc; tic


    % PACKET ERROR RATE BOUNDS - LOWER BOUNDS ON ERROR RATE

    R = 1/2;
    Pev = 10.^([-0.5:-0.25:-1.75,-2:-0.5:-7]);
    uno = ones(size(Pev));

    % normal approximation - expected execution time 15 sec
    SNRdBNA = biawgnPPVbound(n*uno,Pev,R*uno,'normal','error');toc; tic
    % O(n^-2) approximation - expected execution time 120 sec
    SNRdBPPV = biawgnPPVbound(n*uno,Pev,R*uno,'On2','error');toc; tic
    % O(n^-3) approximation - expected execution time 280 sec
    SNRdBPPVb = biawgnPPVbound(n*uno,Pev,R*uno,'On3','error');toc; tic

    save('example0.mat')
else
    load('example0.mat')
end

% display results
close all
figure(1)
set(0,'defaulttextinterpreter','latex')
subplot(1,2,1)
semilogy(EbN0PPV,rhoPPV,'-',EbN0PPVb,rhoPPVb,'--',EbN0NA,rhoNA,'-.')
xlabel('SNR $E_b/N_0$')
ylabel('spectral efficiency $\rho$ [bit/s/Hz]')
title(['n = ' num2str(n) ', Pe = ', num2str(Pe)])
legend({},'interpreter','latex')
legend('$O(n^{-2}\,)$ approx $\quad$','$O(n^{-3}\,)$ approx',...
    'normal approx','Location','Best')
grid
subplot(1,2,2)
semilogy(SNRdBPPV,Pev,'-',SNRdBPPVb,Pev,'--',SNRdBNA,Pev,'-.')
xlabel('SNR $E_b/N_0$')
ylabel('error probability $P_e$')
title(['n = ' num2str(n) ', R = ', num2str(R)])
grid
