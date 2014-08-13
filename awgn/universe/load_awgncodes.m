%
% script that feels out CODES structure for gen_dolinar.m 
%
c = struct('name', 'Hamming', 'ebno', 6.72, 'n', 8, 'k', 4, 'pe', 1e-4, 'fname', 'hamming', 'comment', '');
CODES = c;
c = struct('name', 'Golay', 'ebno', 5.1, 'n', 24, 'k', 12, 'pe', 1e-4, 'fname', 'golay', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Quadratic residue', 'ebno', 4.27, 'n', 48, 'k', 24, 'pe', 1e-4, 'fname', 'qr', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Convolutional (7,1/2)', 'ebno', 4.8, 'n', 48, 'k', 24, 'pe', 1e-4, 'fname', 'conv4', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Convolutional (7,1/2)', 'ebno', 4.72, 'n', 90, 'k', 45, 'pe', 1e-4, 'fname', 'conv3', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Convolutional (7,1/2)', 'ebno', 4.75, 'n', 200, 'k', 100, 'pe', 1e-4, 'fname', 'conv', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Convolutional (7,1/2)', 'ebno', 4.7, 'n', 340, 'k', 170, 'pe', 1e-4, 'fname', 'conv2', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/3', 'ebno', 2.0, 'n', 750, 'k', 250, 'pe', 1e-4, 'fname', 'turbo1', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/3', 'ebno', 1.5, 'n', 1500, 'k', 500, 'pe', 1e-4, 'fname', 'turbo2', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/6', 'ebno', 0.9, 'n', 3000, 'k', 500, 'pe', 1e-4, 'fname', 'turbo3', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/3', 'ebno', 1.0, 'n', 3000, 'k', 1000, 'pe', 1e-4, 'fname', 'turbo4', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/4', 'ebno', 0.8, 'n', 4000, 'k', 1000, 'pe', 1e-4, 'fname', 'turbo5', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/3', 'ebno', 0.85, 'n', 5100, 'k', 1700, 'pe', 1e-4, 'fname', 'turbo7', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/6', 'ebno', 0.58, 'n', 6000, 'k', 1000, 'pe', 1e-4, 'fname', 'turbo6', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Voyager', 'ebno', 2.5, 'n', 8160, 'k', 2*8*223, 'pe', 1e-4, 'fname', 'voyager1', 'comment', 'The exact code: CC(7,1/2)+RS(255,223), Interleaver=2');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/6', 'ebno', 0.3, 'n', 10200, 'k', 1700, 'pe', 1e-4, 'fname', 'turbo8', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/3', 'ebno', 0.6, 'n', 12000, 'k', 4000, 'pe', 1e-4, 'fname', 'turbo9', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Voyager', 'ebno', 2.36, 'n', 16320, 'k', 4*8*223, 'pe', 1e-4, 'fname', 'voyager2', 'comment', 'The exact code: CC(7,1/2)+RS(255,223), Interleaver=4');
CODES = [CODES c];
c = struct('name', 'Galileo HGA', 'ebno', 1.4, 'n', 16320, 'k', 2*8*223, 'pe', 1e-4, 'fname', 'gal_hga', 'comment', 'The exact code: CC(15,1/4)+RS(255,223), Interleaver=2');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/2', 'ebno', 1.03, 'n', 20000, 'k', 10000, 'pe', 1e-4, 'fname', 'turbo11', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Voyager', 'ebno', 2.35, 'n', 20400, 'k', 5*8*223, 'pe', 1e-4, 'fname', 'voyager3', 'comment', 'The exact code: CC(7,1/2)+RS(255,223), Interleaver=5');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/3', 'ebno', 0.4, 'n', 30000, 'k', 10000, 'pe', 1e-4, 'fname', 'turbo12', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/2', 'ebno', 0.8, 'n', 30000, 'k', 15000, 'pe', 1e-4, 'fname', 'turbo12', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Galileo HGA', 'ebno', 1.1, 'n', 32640, 'k', 4*8*223, 'pe', 1e-4, 'fname', 'gal_hga2', 'comment', 'The exact code: CC(15,1/4)+RS(255,223), Interleaver=4');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/4', 'ebno', 0.2, 'n', 40000, 'k', 10000, 'pe', 1e-4, 'fname', 'turbo10', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/3', 'ebno', 0.28, 'n', 45000, 'k', 15000, 'pe', 1e-4, 'fname', 'turbo13', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Turbo R=1/6', 'ebno', -0.1, 'n', 60000, 'k', 10000, 'pe', 1e-4, 'fname', 'turbo14', 'comment', '');
CODES = [CODES c];
c = struct('name', 'Cassini/Pathfinder', 'ebno', 0.9, 'n', 61200, 'k', 5*8*223, 'pe', 1e-4, 'fname', 'cassini', 'comment', 'The exact code: CC(15,1/6)+RS(255,223), Interleaver=5');
CODES = [CODES c];
c = struct('name', 'Galileo LGA', 'ebno', 0.5, 'n', 65280, 'k', 14288, 'pe', 1e-4, 'fname', 'gal_lga', 'comment', 'The exact code: CC(14,1/4)+RSvar(2040,1786)');
CODES = [CODES c];

% These are from T. Richardson
%c = struct('name', 'ME LDPC R=1/2', 'ebno', 2.524, 'n', 640, 'k', 320, 'pe', 1e-4, 'fname', 'ldpc1', 'comment', '');
%CODES = [CODES c];
%c = struct('name', 'ME LDPC R=1/2', 'ebno', 2.021, 'n', 1280, 'k', 640, 'pe', 1e-4, 'fname', 'ldpc2', 'comment', '');
%CODES = [CODES c];
%c = struct('name', 'ME LDPC R=1/2', 'ebno', 1.545, 'n', 2560, 'k', 1280, 'pe', 1e-4, 'fname', 'ldpc3', 'comment', '');
%CODES = [CODES c];
%c = struct('name', 'ME LDPC R=1/2', 'ebno', 1.207, 'n', 5120, 'k', 2560, 'pe', 1e-4, 'fname', 'ldpc4', 'comment', '');
%CODES = [CODES c];
%c = struct('name', 'ME LDPC R=1/2', 'ebno', 0.971, 'n', 10240, 'k', 5120, 'pe', 1e-4, 'fname', 'ldpc5', 'comment', '');
%CODES = [CODES c];

% Fig. 3 of Lee-O'Sullivan "Alg. soft decision decoding of Hermitian Codes",TIT-56, 6, June 2010
% The code is BPSK modulated AG-code [64,32]_{4^2} based on a Hermitian curve: y^4+y-x^5 over F_{4^2}
% 	genus = 6, Line bundle = L(uP_\infty), u=37; mindist = n-u = 27;
%
% Note: the resulting code efficiency is 
c = struct('name', 'Hermitian curve [64,32] (SDD)', 'ebno', 5.877, 'n', 4*64, 'k', 4*32, 'pe', 1e-4, 'fname', 'herm1', 'comment', '');
CODES = [CODES c];


% Fig. 2 and 3 of 
% RS[31,27,5]
c = struct('name', 'Reed-Solomon (SDD)', 'ebno', 5.24, 'n', 31*8, 'k', 27*8, 'pe', 1e-4, ...
	'fname', 'rs3', 'comment', 'ML decoder');
CODES = [CODES c];
% RS[31,25,7]
c = struct('name', 'Reed-Solomon (SDD)', 'ebno', 5.57, 'n', 31*8, 'k', 25*8, 'pe', 1e-4, ...
	'fname', 'rs3', 'comment', 'SUB decoder');
CODES = [CODES c];
% RS[63,59,5]
c = struct('name', 'Reed-Solomon (SDD)', 'ebno', 6.34, 'n', 63*8, 'k', 59*5, 'pe', 1e-4, ...
	'fname', 'rs3', 'comment', 'SUB decoder');
CODES = [CODES c];
% RS[127,125,3]
c = struct('name', 'Reed-Solomon (SDD)', 'ebno', 7.56, 'n', 127*8, 'k', 125*8, 'pe', 1e-4, ...
	'fname', 'rs4', 'comment', 'ML decoder');
CODES = [CODES c];

% Fig. 4 of Gross,Kschischiang, Koetter, Gulak "Simul. Res. for Alg. Soft Dec. Dec. of RS Codes"
% RS code [255, 239] over BPSK modulated AWGN; ebno = 6.675; pe = 1e-4.
c = struct('name', 'Reed-Solomon (SDD)', 'ebno', 6.675, 'n', 8*255, 'k', 8*239, 'pe', 1e-4, ...
	'fname', 'rs1', 'comment', 'Koetter-Vardy decoder');
CODES = [CODES c];

% Fig. 2 of Koetter-Vardy patent appl. RS(256, 144, 113) concatenated with 8-to-9 binary code (BPSK);
% SNR = 3.85 dB => Ebno = 3.85 db (since overall rate is 1/2)
c = struct('name', 'Reed-Solomon (SDD)', 'ebno', 3.85, 'n', 144*8*2, 'k', 144*8, 'pe', 1e-4, ...
	'fname', 'rs2', 'comment', 'Koetter-Vardy decoder');
CODES = [CODES c];



% Also from Koetter-Vardy patent appl. BPSK mod. BCH code [127, 71, 19] + KV alg.
c = struct('name', 'BCH (Koetter-Vardy)', 'ebno', 4.571, 'n', 127, 'k', 71, 'pe', 1e-4, ...
	'fname', 'bch1', 'comment', '');
CODES = [CODES c];

%
%
% Attention: I dropped this point (it has norm. rate \approx 1.03)
%c = struct('name', 'BCH (ML decoding)', 'ebno', 3.186, 'n', 127, 'k', 71, 'pe', 1e-4, ...
%	'fname', 'bch2', 'comment', '');
%CODES = [CODES c];


%
% Polar code + CRC + List decoder (Vardy-Tal)
%
% Note: the k is the actual number of information bits (thus the actual polar code operates with k = 1024+16 bits of CRC)
% 
%c = struct('name', 'Polar+CRC R=1/2 (List dec.)', 'ebno', 1.519, 'n', 2048, 'k', 1024, 'pe', 1e-3, 'fname', 'polar1', 'comment', 'The exact code: Vardy-Tal list-decoded polar+CRC');
c = struct('name', 'Polar+CRC R=1/2 (List dec.)', 'ebno', 1.8283, 'n', 2048, 'k', 1024, 'pe', 1e-4, 'fname', 'polar1', 'comment', 'The exact code: Vardy-Tal list-decoded polar+CRC');
CODES = [CODES c];

c = struct('name', 'Polar+CRC R=1/2 (List dec.)', 'ebno', 2.09063, 'n', 1024, 'k', 512, 'pe', 1e-4, 'fname', 'polar2', 'comment', 'The exact code: Vardy-Tal list-decoded polar+CRC');
CODES = [CODES c];



% These are from T. Richardson's program via scaling
if(universe_only)
	addpath ./me_ldpc
	for n=[500:100:1000 1200:200:2000 2400:400:4000 5000:1000:10000 10000:2500:20000];
		[ss ee] = me_snr(n, 1/2, 1e-4);
		ebno = 10*log10(ee);
		c = struct('name', 'ME LDPC R=1/2 (BP)', 'ebno', ebno, 'n', n, 'k', n/2, ...
			'pe', 1e-4, 'fname', '...', 'comment', '');
		CODES = [CODES c];
	end
	rmpath ./me_ldpc
end
