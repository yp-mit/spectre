% see comments in plot_universe.m

% Universe only assumed
universe_only = 1;
skip_cc48 = 1;

load_awgncodes

% These are to sort codes in classes
%          b     blue          .     point              -     solid
%          g     green         o     circle             :     dotted
%          r     red           x     x-mark             -.    dashdot 
%          c     cyan          +     plus               --    dashed   
%          m     magenta       *     star             (none)  no line
%          y     yellow        s     square
%          k     black         d     diamond
%                              v     triangle (down)
%                              ^     triangle (up)
%                              <     triangle (left)
%                              >     triangle (right)
%                              p     pentagram
%                              h     hexagram

CLASSES = struct('class', [], ...
	'plot', {'-kx', '-b+', '-g*', '-rs', '-kd', '-b^', '-gp', '-rh', '--ko', '--r+', '-bo', '-go', '-ro', '-ko',...
		 '-rx', '-bx', '-gx', '-ks', '-bs', '-gs'}, ...
	'idxs', [0]);

addpath ./../

N_codes = size(CODES,2);
figure(1);

last_class = 0;

for idx = 1:N_codes;
	Nmin = 10;
	Nmax = max(500, ceil(1.1*CODES(idx).n));
	epsil = CODES(idx).pe;
	code_rate = (CODES(idx).k)/(CODES(idx).n);
	ebno = 10^((CODES(idx).ebno)/10);
	P = 2*code_rate*ebno;
	C = cap_awgn(P);
	Ns_cap = floor(linspace(0, Nmax, 9));
	Capr = C + 0 .* Ns_cap;
	
	Ns_norm = floor(linspace(80, Nmax));
	nrm = normapx_awgn(Ns_norm, epsil, P);
	Ns_norm_add = floor(linspace(10,80));
	nrm_add= normapx_awgn(Ns_norm_add, epsil, P);

	found = 0;

	for cc=1:last_class;
		if (strcmp(CLASSES(cc).class, CODES(idx).name));
			CLASSES(cc).idxs = [CLASSES(cc).idxs idx];
			found = 1;
			break;
		end
	end
	
	if (~found)
		% Exclude a few elements
		if (strcmp(CODES(idx).name, 'Hamming')) || ...
			(strcmp(CODES(idx).name, 'Golay')) || ...
			(strcmp(CODES(idx).name, 'Quadratic residue')) || ...
			((strcmp(CODES(idx).name, 'Convolutional (7,1/2)'))) ||...
			((skip_cc48) && strcmp(CODES(idx).name, 'Convolutional (7,1/2)') && (CODES(idx).n == 48))
			disp(sprintf('-- Summary: skipping: %s (%d, %d).', ...
					CODES(idx).name, CODES(idx).n, CODES(idx).k));
		else
			last_class = last_class + 1;
			if(last_class > size(CLASSES, 2))
				disp('ERROR: too many classes!');
				error('plot_universe');
			end
			CLASSES(last_class).class = CODES(idx).name;
			CLASSES(last_class).idxs = idx;
		end
	end
end

% Generate a summary figure
base_pe = 1e-4;

figure(1); fig1=gcf; axes('FontSize', 14);
leg = {''};
for cc=1:last_class;
	cls_size = size(CLASSES(cc).idxs, 2);
	blocklens=[];
	dbgaps = [];
	for kk=1:cls_size;
		idx = CLASSES(cc).idxs(kk);
		epsil = CODES(idx).pe;
		code_rate = (CODES(idx).k)/(CODES(idx).n);
		ebno = 10^((CODES(idx).ebno)/10);
		P = 2*code_rate*ebno;
		if(epsil ~= base_pe)
			disp(sprintf('ERROR: epsil != base_pe for the code with idx=%d\n', idx));
			error('plot_universe');
		end
		bllen = CODES(idx).n;
		Popt = optpower(epsil, bllen, code_rate);
		dbgap = 10*log10(P) - 10*log10(Popt);
		blocklens = [blocklens bllen];
		dbgaps = [dbgaps dbgap];
	end
	plotstr = CLASSES(cc).plot;
	if(cls_size == 1)
		if (strcmp(plotstr(1:2), '--'))
			plotstr = plotstr(3:end);
		elseif (plotstr(1)=='-')
			plotstr = plotstr(2:end);
		end
	end
	semilogx(blocklens, dbgaps, plotstr, 'LineWidth', 1.0, 'MarkerSize', 7.0);
	if(cc == 1)
		hold on;
	end
	leg(cc) = {CLASSES(cc).class};
end
%ylim([0.5 1.1]); 
grid on;
title(sprintf('dB gaps to finite blocklength fund. limit for codes over AWGN, Pe=%g', base_pe));
xlabel('Blocklength, n'); ylabel('dB gap');
legend(leg); legend('Location', 'SouthWest');
set(fig1, 'PaperPositionMode', 'manual');
set(fig1, 'PaperPosition', [0 0 12 9]);
figure(fig1); 
	xlim([10^2 10^5]);
	legend('Location', 'EastOutside');
	print('-depsc2', 'output/universe_db.eps');


