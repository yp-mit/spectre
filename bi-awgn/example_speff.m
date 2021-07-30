% SPECTRAL EFFICIENCY CURVES - UPPER BOUNDS ON RATE

if(0) % set to (1) if you wish to run the converse_mc function

    n = 1e3;
    Pe = 1e-5;
    SNRdB = -15:7;
    uno = ones(size(SNRdB));tic

    % normal approximation - expected execution time 0.5 sec
    tic; rhoNA = 2*converse_mc(n*uno,Pe*uno,SNRdB,'normal'); % spectral efficiency 
    EbN0NA = SNRdB-10*log10(rhoNA); toc
    % O(n^-2) approximation - expected execution time 10 sec
    tic; rhoPPV2 = 2*converse_mc(n*uno,Pe*uno,SNRdB,'On2'); % spectral efficiency 
    EbN0PPV2 = SNRdB-10*log10(rhoPPV2); toc
    % O(n^-3) approximation - expected execution time 20 sec
    tic; rhoPPV3 = 2*converse_mc(n*uno,Pe*uno,SNRdB,'On3'); % spectral efficiency 
    EbN0PPV3 = SNRdB-10*log10(rhoPPV3); toc
    % full approximation - expected execution time 10 min
    tic; rhoPPV = 2*converse_mc(n*uno,Pe*uno,SNRdB,'full'); % spectral efficiency 
    EbN0PPV = SNRdB-10*log10(rhoPPV); toc

    % save('example_speff.mat')
else
    load('example_speff.mat')
end

% display results
close all
figure(1)
set(0,'defaulttextinterpreter','latex')
semilogy(EbN0PPV3,rhoPPV3,'--', EbN0PPV,rhoPPV,'-',EbN0PPV2,rhoPPV2,'--',EbN0NA,rhoNA,'-.x')
xlabel('SNR $E_b/N_0$')
ylabel('spectral efficiency $\rho$ [bit/s/Hz]')
title(['n = ' num2str(n) ', Pe = ', num2str(Pe)])
legend({},'interpreter','latex')
legend('$O(n^{-3}\,)$ approx $\quad$','full approx','$O(n^{-2}\,)$ approx',...
    'normal approx','Location','Best')
grid
