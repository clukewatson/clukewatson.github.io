function [si_j] = ind_share(argMU,argPSI)
%{
%%%THIS IS FOR RCLOGIT!
THETANL	= thetaNL_init;
MU 	= mu_init;
PSI 	= psifunc(THETANL);
%}

% *********************************************************************
% This uses the Overflow Safe Method, as advocated by pyblp
% *********************************************************************

% Utilities
V 	= argMU + argPSI;

% Perform the Safety
mxV 	= max(V);
Vp 	= V - mxV;
scale 	= exp(-mxV);

% Numerator
num 	= exp(Vp);

% Denomenator
den 	= ( scale + sum(num,1) );

% *********************************************************************
% Calculate Share
% *********************************************************************
% NOTE: set NaN values to zero 
%	this is a little hacky

si_j 	= num ./ (den);
si_j(isnan(si_j)) = 0;


end