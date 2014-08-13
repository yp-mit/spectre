function [snr ebno] = me_snr(n, R, Pe)
% Gets the snr (and ebno) for the Multi-Edge LDPC of rate R and blocklength n
% note: ebno is defined similar to the gen_dolinar:
% SNR = 2*ebno*R   <=> ebno = SNR/2/R
%
% Note: all output values are NOT in db

% hardcoded constraints in compute_ab
k = n*R;
if (k < 250) || (k > 10000) || (R < 25/74) || (R>25/32)
	disp(sprintf('me_snr: n=%d, R=%g, Pe=%g: ERROR: (k,R) out of range', n, R, Pe));
	snr = NaN;
	ebno = NaN;
	return;
end
[a b] = compute_ab(k, R);
snr_db = (-norminv(Pe))/a + b;

snr = 10^(snr_db/10);
ebno = snr/2/R;

