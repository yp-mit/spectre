% Script to show various bounds and normal approximations for the energy-per-bit

% First plot is for the AWGN channel

epsil = 1e-3;
Es = 10.^linspace(-.3,6);       
Lms_conv = energy_awgn_conv(Es, epsil);
Lms_na = max(energy_awgn_normapx(Es, epsil), 0);
Es_na = Es(Lms_na > 10); Lms_na = Lms_na(Lms_na > 10);

Es_opt = [min(Es) max(Es)];
Lms_opt = Es_opt./log(2);

Lms_ach = 10.^linspace(0, 4);
Es_ach = energy_awgn_ach(Lms_ach, epsil);

col_red = [1 .2 0];
col_green = [0 .6 0];
col_blue = [0 .1 1];

figure; fig1=gcf; axes('FontSize', 14);
semilogx(Lms_ach, 10*log10(Es_ach./Lms_ach), 'LineWidth', 1.0, 'Color', col_blue); hold on;
semilogx(Lms_na, 10*log10(Es_na./Lms_na), 'k--', 'LineWidth', 1.0); hold on;
semilogx(Lms_conv, 10*log10(Es./Lms_conv), 'LineWidth', 1.0, 'Color', col_red); hold on;
semilogx(Lms_opt, 10*log10(Es_opt./Lms_opt), '--', 'LineWidth', 1.0, 'Color', col_red); hold on;

xlabel('Information bits, k');
ylabel('Eb/No, dB');
title(sprintf('Energy per bit vs. length of the message (AWGN channel); \\epsilon=%g', epsil));
legend('Achievability', 'Normal approximation', ...
			'Converse', 'Shannon limit, -1.59 dB');

xlim([1 10^6]);

%set(fig1, 'PaperPositionMode', 'manual');
%set(fig1, 'PaperPosition', [0 0 12 9]);
%figure(fig1);
%print -depsc2 energy_awgn.eps


figure; fig2=gcf; axes('FontSize', 14);

EE = 10.^linspace(2,5,30);
Lm_na = max(energy_nocsi_normapx(EE, epsil),0);
Lm_ach_ht = energy_nocsi_ach_ht(EE,epsil);

disp('Starting computation of 30 points of converse. This may take several hours...');

%Lm_conv = energy_nocsi_conv(EE,epsil);
Lm_conv = 0*EE;

EE_csir_na = 10.^linspace(1.5,5,100);
Lm_csir_na = max(energy_awgn_normapx(EE_csir_na, epsil), 0);

Lms_opt = [1 2e5];
Es_opt = Lms_opt * log(2);

semilogx(Lm_ach_ht, 10*log10(EE./Lm_ach_ht), 'LineWidth', 1.0, 'Color', col_blue); hold on;
semilogx(Lm_na, 10*log10(EE./Lm_na), 'k-', 'LineWidth', 1.0); hold on;
semilogx(Lm_conv, 10*log10(EE./Lm_conv), 'LineWidth', 1.0, 'Color', col_red); hold on;
semilogx(Lm_csir_na, 10*log10(EE_csir_na./Lm_csir_na), 'k--','LineWidth', 1.0); hold on;
semilogx(Lms_opt, 10*log10(Es_opt./Lms_opt), 'r--', 'LineWidth', 1.0); hold on;
xlabel('Information bits, k');
ylabel('Eb/No, dB');
title(sprintf('Energy per bit vs. length of the message (Rayleigh fading); \\epsilon=%g', epsil));
legend('noCSI: Achievability (ach\_ht)', 'noCSI: Normal approximation', ...
			'noCSI: Converse', 'CSIR: Normal approximation', 'Shannon limit, -1.59 dB');
xlim([1 2e5]);
ylim([-2 10]);
%set(fig2, 'PaperPositionMode', 'manual');
%set(fig2, 'PaperPosition', [0 0 12 9]);
%figure(fig2);
%print -depsc2 energy_nocsi.eps
