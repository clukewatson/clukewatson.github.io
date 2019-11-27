function s_j = mrktshare(argMU,argPSI)
%{
THETA = theta_init;
MU = mu_init;
PSI = psifunc(THETA);
%}

s_j = mean(ind_share(argMU,argPSI),2);

end