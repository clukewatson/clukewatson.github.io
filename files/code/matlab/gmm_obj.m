function [m] = gmm_obj(BETA)
%{
BETA = B_ols
%}
global J W y X Z;

resid = y - X*BETA;
ZR    = (resid'*Z);
m     = ZR*W*ZR'/J;
%m =  ((y - X*BETA)'*Z)*W*((y - X*BETA)'*Z)'/J;

end