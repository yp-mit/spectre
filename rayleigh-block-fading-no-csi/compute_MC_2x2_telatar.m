function R=compute_MC_2x2_telatar(snrdB,T,L,epsilon)

Mt=2;
Mr=2;

prec=round(log2(100/epsilon));

numIt=L+1;
pow_all=zeros(1,L);
Ralt1=zeros(numIt,1);
Ralt2=zeros(numIt,1);


for ii=1:1:numIt, % try all points according to generalized Telatar conjecture
    numIt-ii
    bstring=dec2bin(2^(ii-1)-1,L);
    
    for k=1:1:L,
        pow_all(k)=str2num(bstring(k))*0.5;        
    end
    
    filename = ''; %create filename of raw data
    
    Ralt1(ii)=MC_USTM_2x2(snrdB,T,L,1,epsilon,prec,pow_all,filename)
    Ralt2(ii)=MC_USTM_2x2(snrdB,T,L,2,epsilon,prec,pow_all,filename)
    
end

Ralt1_max=max(Ralt1);
Ralt2_max=max(Ralt2);

R=min(Ralt1_max,Ralt2_max);

