%% 1) SPECTRAL EFFICIENCY CURVES
disp(' ')
disp('SPECTRAL EFFICIENCY CURVES')

Pe = 1e-5;
data = {... % n, om_NA, om_PPV
    {{2e2},{-5:10},{-7:7}},...
    {{1e3},{-12:10},{-15:9}},...
    {{1e4},{-22:10},{-25:10}},...
    {{1e5},{-32:10},{-30:10}},...
    {{1e6},{-32:10},{-32:10}}...
};
snr_dB = -32:15;


if (0) % set to 1 if you want to generate values
    
    % NA bound 
    % expected execution time 8 sec
    rho_NA = nan(length(data),length(snr_dB));
    EbN0_NA = nan(length(data),length(snr_dB));  tic
    for k=1:length(data)
        n = data{k}{1}{1};
        Om_dB = data{k}{2}{1};
        uno = ones(size(Om_dB));
        R_NA = biawgnPPVbound(n*uno,Pe*uno,Om_dB,'normal'); 
        rho_NA(k,ismember(snr_dB,Om_dB)) = 2*R_NA;
        EbN0_NA(k,ismember(snr_dB,Om_dB)) = Om_dB-10*log10(2*R_NA);
    end; toc

    % PPV bound 
    % expected execution time 250 sec (4 minutes)
    rho_PPV = nan(length(data),length(snr_dB));
    EbN0_PPV = nan(length(data),length(snr_dB)); tic
    for k=1:length(data)
        n = data{k}{1}{1};
        Om_dB = data{k}{3}{1};
        uno = ones(size(Om_dB));
        R_PPV = biawgnPPVbound(n*uno,Pe*uno,Om_dB,'On2'); 
        rho_PPV(k,ismember(snr_dB,Om_dB)) = 2*R_PPV;
        EbN0_PPV(k,ismember(snr_dB,Om_dB)) = Om_dB-10*log10(2*R_PPV);
    end; toc

    save('example1.mat')
else
    load('example1.mat')
end

% Figures

figure(1)
set(0,'defaulttextinterpreter','latex')
legend({},'interpreter','latex')
subplot(1,1,1)
semilogy(EbN0_NA.',rho_NA.','-.')
hold on
semilogy(EbN0_PPV.',rho_PPV.')
hold off
xlabel('SNR $E_b/N_0$')
ylabel('spectral efficiency $\rho$ [bit/s/Hz]')
grid
