function djac = jacob_rclogit(argMU,argPSI)
%{
THETANL = theta_init;
MUold 	= mu_init;
PSIold 	= psifunc(THETANL);
%}

global P yi vi Xc; 

% Get Trial Individual Level Shares 
[si_j] 	= ind_share(argMU,argPSI);

% *********************************************************************
% Calculate GMM Gradient Parts
% *********************************************************************

d = cell(1,7);
for i = 1 : 7

% Derivative wrt Alpha -- coef on price
  if i == 1; d{i}=DSDP(si_j,P,yi) ; end

% Derivatives wrt Sigma -- coef on X
  if i == 2; d{i}=DSDX_VI(si_j,vi(:,1),Xc(:,1)); end
  if i == 3; d{i}=DSDX_VI(si_j,vi(:,2),Xc(:,2)); end
  if i == 4; d{i}=DSDX_VI(si_j,vi(:,3),Xc(:,3)); end
  if i == 5; d{i}=DSDX_VI(si_j,vi(:,4),Xc(:,4)); end
  if i == 6; d{i}=DSDX_VI(si_j,vi(:,5),Xc(:,5)); end
   
% Derivative wrt MU -- mean value   
  if i == 7; d{i}=DSDX_MU(si_j); end  

end

% Jacobian
djac = (-1) * ([d{7}] \ [d{1:6}]);

end
