function dsdx_v = DSDX_VI(SJ,VIC,XC)
%{
This program calculates ds/dsigma.
MU 	= mu_init;
PSIold 	= psifunc(thetaNL_init);
THETANL = thetaNL_init;
YI 	= yi;
SJ 	= SH_j;
PP 	= P;
VIC 	= vi(:,1);
XC 	= Xc(:,1);
%}

%%%% Derivative wrt sigma_v_l
%
dsdx_v = mean( VIC'.*SJ.*( XC - ( sum(XC.*SJ,1) ) ), 2);

end

