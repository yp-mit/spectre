function plot_all(distortion, excess_prob, N);

path(path, '../../lib');

clear global;
if nargin<1
	distortion = 0.1;
	excess_prob = 0.01;
	N=10:5:500;
end

for k=1:length(N); 
	out_rogers(k) = Rstar(distortion, N(k), excess_prob, 'rogersa');
	out_conv(k) = Rstar(distortion, N(k), excess_prob, 'spherecoveringc');
	out_norm(k) = Rstar(distortion, N(k), excess_prob, 'normal');
end;

figure(1); plot(N, out_rogers, 'b', N, out_norm, 'k', N, out_conv, 'r'); 
xlabel('Blocklength, n'); ylabel('Compression rate, bit / source sample');
title(sprintf('Lossy compression of iid standard normal source. Distortion=%g at excess=%g%%', distortion, 100*excess_prob));

legend('Rogers (achievability)', 'Normal approximation', 'Sphere covering + HT (converse)'); 
