
function R=outage_2x2_telatar(tfdiv,prec,epsilon,snrdb)
  %
  % optimize outage mi over all input covariance matrices of the form given in [3, eq. ()] (generalized Telatar conjecture)
  % tfdiv: number of time-frequency diversity branches
  % prec: log2 of the number of samples used in the Monte-Carlo step
  % epsilon: outage probability
  % snrdb: snr in db.
  
mt=2;
mr=2;

numIt=tfdiv+1;
pow_all=zeros(1,tfdiv);
R=zeros(numIt,1);

for ii=[numIt],
    numIt-ii;
    bstring=dec2bin(2^(ii-1)-1,tfdiv);
    
    for k=1:1:tfdiv,
        pow_all(k)=str2num(bstring(k))*0.5;        
    end
    
    tic
    R(ii)=outage_mi(snrdb,epsilon,mt,mr,tfdiv,pow_all,prec)
    elapsed_time=toc/60;
    msg=sprintf('time needed to generate samples: %f min', elapsed_time);
    disp(msg); 
end
