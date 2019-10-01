clear all
close all
close hidden
warning off all
clc

% Program Parameters
cd 'mypath';

%**************************************************************************
% GOAL:
% Replicate stata commands:
%   reg y x* w
%   ivreg2 gmm y x* (w = z*)
% 
% Follow the same DGP and you should be quite close.
%**************************************************************************

%**************************************************************************
%Define globals
%**************************************************************************
global N J W y X Z;
rng(52)

%**************************************************************************
% Parameters
%**************************************************************************
N = 1;
J = 10000;
B0 = [ 1 3 ]

mytolx        = 1e-6;
mytolfun      = 1e-6;
mymaxiters    = 5*10^5;
mymaxfunevals = 4000;

%**************************************************************************
% Random Numbers
%**************************************************************************

mu 	= [0 0];
sigma 	= [1 0.5; 0.5 1];
R 	= mvnrnd(mu,sigma,J);
e_w 	= R(:,1);
e_y 	= R(:,2);
% plot(R(:,1),R(:,2),'+')
% cov(e_w,e_y)

Rz 	= mvnrnd(mu,sigma,J);
z3 	= Rz(:,1);
e_u3 	= Rz(:,2);

z1 	= randn(J,N);
z2 	= randn(J,N);
w 	= 1 + z1 + z2 + z3 + e_w;
y 	= B0(1)*1 + B0(2)*w + e_y + e_u3;
X 	= [ones(J,1) w];
Xexog 	= [ones(J,1)];
Z 	= [Xexog z1 z2];
Zb 	= [Xexog z1 z2 z1.*z2 z1.*z1 z2.*z2];

% True
B0

% OLS
B_ols 	= ((X'*X)\(X'*y))

% IV regressions
PZ 	= Z*inv(Z'*Z)*Z';
B_2sls 	= ((X'*PZ*X) \ (X'*PZ*y))
% == ((X'*(Z*inv(Z'*Z)*Z')*X) \ (X'*(Z*inv(Z'*Z)*Z')*y))'
% other stuff
% FS 		= Z*((Z'*Z) \ Z'*w);
% XX 		= [Xexog FS];
% B_2sls_v2 	= ((XX'*XX)\XX'*y)

% GMMIV Two Step using Identity Weight Matrix
W 	= eye(size(Z,2));
Bi_gmm1 = ((X'*Z*W*Z'*X) \ (X'*Z*W*Z'*y))
ZX 	= (Z'*X)'/J;
r 	= y - X*Bi_gmm1;
Zu 	= (Z.*(r*ones(1,size(Z,2))));
ZuuZ 	= Zu'*Zu/J;
W 	= inv(ZuuZ);
Bi_gmm2 = ((X'*Z*W*Z'*X) \ (X'*Z*W*Z'*y))


% GMMIV Two Step using 2SLS Initial Weight
W 	= inv(Z'*Z/J);
Bz_gmm1 = ((X'*Z*W*Z'*X) \ (X'*Z*W*Z'*y))
ZX 	= (Z'*X)'/J;
r 	= y - X*Bz_gmm1;
Zu 	= (Z.*(r*ones(1,size(Z,2))));
ZuuZ 	= Zu'*Zu/J;
W 	= inv(ZuuZ);
Bz_gmm2 = ((X'*Z*W*Z'*X) \ (X'*Z*W*Z'*y))



% % Digression
% %   Skipping the Js does not change anything
% %       for estimating the Beta 
% %	  since they cancel out
% % 
% W = inv(Z'*Z);
% B_gmm1_noJ = ((X'*Z*W*Z'*X) \ (X'*Z*W*Z'*y));
% ZX = (Z'*X)';
% r = y - X* B_gmm1_noJ;
% Zu = (Z.*(r*ones(1,size(Z,2))));
% ZuuZ = Zu'*Zu;
% W = inv(ZuuZ);
% B_gmm2_noJ = ((X'*Z*W*Z'*X) \ (X'*Z*W*Z'*y))


% GMM NonLinear Search
B_init = B_ols;
options=optimset('TolFun',mytolfun, ...
        'TolX',mytolx, ...
        'Display','iter', ...
        'MaxIter',mymaxiters, ...
        'MaxFunEvals',mymaxfunevals);

% Quasi-Newton Method (Gradient Based for smooth Obj Functions)    
%   using finite differences (all defaults)
[B_gmmNL_QN] = fminunc('gmm_obj',B_init,options);

% Non-derivative method for any obj function
%[B_gmmNL_NM] = fminsearch('gmm_obj',B_init,options);

% % As Annonymous Functions
% gmm_obj_annon = @(BETA) (((y - X*BETA)'*Z)*W*((y - X*BETA)'*Z)'/J);
% [B_gmmNL_QN_an] = fminunc(gmm_obj_annon,B_init,options);
% [B_gmmNL_NM_an] = fminsearch(gmm_obj_annon,B_init,options);
% disp([B_gmmNL_QN B_gmmNL_NM B_gmmNL_QN_an B_gmmNL_NM_an])


% Iterated GMM (Code Based on Hanson & Lee 2019)
tolerance = 1e-9;
maxit = 1e+3;

b = Bz_gmm2;
for iter = 1:maxit
   e = y - X*b;
   w = zeros(l,l);
   ze = Z.*repmat(e,1,l);
   w = (ze'*ze)/n;   
   b = (zx'/w*zx)\(zx'/w*zy);
   db = b - b1;
   if norm(db) < tolerance
       break
   end   
   b1 = b;
   
   if iter == maxit
       b = NaN;
       s = NaN;
       V = NaN;
       sw = NaN;
       Vw = NaN;
       s0 = NaN;
       V0 = NaN;
       J = NaN;
       pv = NaN;
       return
   end
end
B_igmm = b;


disp([B0 B_ols B_iv B_2sls Bi_gmm1 Bz_gmm1 B_gmmNL_QN B_igmm])



% OLS Standard Errors
r_ols = y - X*B_ols;
V_ols = (r'*r / J)*inv(X'*X);
se_ols = sqrt(diag(V_ols));
t_ols = B_ols ./ se_ols;
ols_out = [B_ols, se_ols, t_ols]

% Robust GMM Standard Errors
resid = y - X* B_igmm;
ZX = (Z'*X)'/J;
DELTA = (ZX*W*ZX') \ (ZX*W);
Zu    = (Z.*(resid*ones(1,size(Z,2))));
ZuuZ = Zu'*Zu/J;
V_gmm = (DELTA*ZuuZ*DELTA')/J;
se_gmm = full(sqrt(diag(V_gmm)));
t_gmm = B_igmm ./ se_gmm;
gmm_out = [B_igmm, se_gmm, t_gmm]

% % Digression -- Js still don't matter!
% r = y - X*B_gmm0;
% Zu = (Z.*(r*ones(1,size(Z,2))));
% ZuuZ = Zu'*Zu;
% W = inv(ZuuZ);
% ZX = (Z'*X)';
% DELTA = (ZX*W*ZX') \ (ZX*W);
% Zu    = (Z.*(resid*ones(1,size(Z,2))));
% ZuuZ = Zu'*Zu;
% V_gmm = (DELTA*ZuuZ*DELTA');
% se_gmm = full(sqrt(diag(V_gmm)));
% t_gmm = B_gmm1 ./ se_gmm;
% gmm_out2 = [B_gmm1, se_gmm, t_gmm]


