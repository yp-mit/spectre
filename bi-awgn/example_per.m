% PACKET ERROR RATE BOUNDS - LOWER BOUNDS ON ERROR RATE

if(0) % set to (1) if you wish to run the converse_mc function

    n = 1e3;
    R = 1/2;
    Pev = 10.^([-0.5:-0.25:-1.75,-2:-0.5:-7]);
    uno = ones(size(Pev));

    % normal approximation - expected execution time 15 sec
    SNRdBNA = converse_mc(n*uno,Pev,R*uno,'normal','error');toc; tic
    % O(n^-2) approximation - expected execution time 120 sec
    SNRdBPPV = converse_mc(n*uno,Pev,R*uno,'On2','error');toc; tic
    % O(n^-3) approximation - expected execution time 280 sec
    SNRdBPPVb = converse_mc(n*uno,Pev,R*uno,'On3','error');toc; tic

    save('example_per.mat')
else
    load('example_per.mat')
end

% display results
close all
figure(1)
set(0,'defaulttextinterpreter','latex')
semilogy(SNRdBPPV,Pev,'-',SNRdBPPVb,Pev,'--',SNRdBNA,Pev,'-.')
xlabel('SNR $E_b/N_0$')
ylabel('error probability $P_e$')
title(['n = ' num2str(n) ', R = ', num2str(R)])
legend({},'interpreter','latex')
legend('$O(n^{-2}\,)$ approx $\quad$','$O(n^{-3}\,)$ approx',...
    'normal approx','Location','Best')
grid
