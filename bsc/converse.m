function Lms = converse(Ns, delta, epsil)

shutup = 0;
if max(size(Ns,2)) > 1
	shutup = 1;
end

Lms = zeros(size(Ns));

for idx=1:length(Ns)
	n = Ns(idx);
	conv = betanq(epsil, n, delta, shutup);
	Lms(idx) = -conv;
end
