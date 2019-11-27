function ds_dmu = DSDX_MU(SJ)
%{
This program calculates ds/dmu.
MU 	= mu_init;
PSI 	= psifunc(thetaNL_init);
THETANL = thetaNL_init;
YI 	= yi;
[SJ] 	= ind_share(MU,PSI,THETANL);
PP 	= P;
%}

% Helpers
N = size(SJ,2);
M = ( SJ * SJ' ) / N ;
D = mean(SJ,2) ;

%%%% Derivative wrt mu
%
ds_dmu = diag(D) - M ;

end


