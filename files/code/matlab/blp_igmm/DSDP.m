function dsdp = DSDP(SJ,PP,YI)
%{
This program calculates ds/dalpha.
THETANL = thetaNL_init;
MU 	= mu_init;
PSIold 	= psifunc(THETANL);
YI 	= yi;
[SJ] 	= ind_share(MU,PSI,THETANL);
PP 	= P;
%}

%%%% Derivative wrt alpha
%
dsdp = mean( (1./YI)'.*SJ.*( PP - (sum(PP.*SJ,1)) ) ,2);

end


