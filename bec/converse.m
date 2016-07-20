function Lms = converse(Ns, delta, epsil);
% This is just a vectorization of the the converse_spec2.m. See comments there

if(max(size(Ns))>1)
	shutup=1;
else
	shutup=0;
end

Lms = [];

for n = Ns;
	Lms = [Lms converse_spec2(n, delta, epsil, [], [], shutup)];
end
