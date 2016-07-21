% SPECTRAL EFFICIENCY CURVES - UPPER BOUNDS ON RATE

if(0) % set to (1) if you wish to run the converse_mc function

    n = 1e3;
    Pe = 1e-5;
    SNRdB = -15:10;
    uno = ones(size(SNRdB));tic

    % normal approximation - expected execution time 2 sec
    rhoNA = 2*converse_mc(n*uno,Pe*uno,SNRdB,'normal'); % spectral efficiency 
    EbN0NA = SNRdB-10*log10(rhoNA);toc; tic
    % O(n^-2) approximation - expected execution time 90 sec
    rhoPPV = 2*converse_mc(n*uno,Pe*uno,SNRdB,'On2'); % spectral efficiency 
    EbN0PPV = SNRdB-10*log10(rhoPPV);toc; tic
    % O(n^-3) approximation - expected execution time 120 sec
    rhoPPVb = 2*converse_mc(n*uno,Pe*uno,SNRdB,'On3'); % spectral efficiency 
    EbN0PPVb = SNRdB-10*log10(rhoPPVb);toc; tic

    save('example_speff.mat')
else
    load('example_speff.mat')
end

% display results
close all
figure(1)
set(0,'defaulttextinterpreter','latex')
semilogy(EbN0PPV,rhoPPV,'-',EbN0PPVb,rhoPPVb,'--',EbN0NA,rhoNA,'-.')
xlabel('SNR $E_b/N_0$')
ylabel('spectral efficiency $\rho$ [bit/s/Hz]')
title(['n = ' num2str(n) ', Pe = ', num2str(Pe)])
legend({},'interpreter','latex')
legend('$O(n^{-2}\,)$ approx $\quad$','$O(n^{-3}\,)$ approx',...
    'normal approx','Location','Best')
grid
