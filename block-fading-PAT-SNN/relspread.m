function rs = relspread(data)
%RELSPREAD returns the relative spread of the data expressed in terms of
%interquartile range over the median. 

q = quantile(data,[.25 .50 .75]);   % the quartiles of eps_trial
iqr_eps = q(3) - q(1);              % interquartile range (measure of spread)
med_eps = q(2);                     % median (measure of location)
rs = iqr_eps / med_eps;             % relative spread