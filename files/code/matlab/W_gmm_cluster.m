function [b,s,V,sw,Vw,s0,V0,J,pv,tw,t] = W_gmm_cluster(y,x,ww,z,mem)
%{
y = y;
ww = W;
x = X;
z = ZT;
mem = Grp_Idx;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated on Sept 22, 2019 by CLW
% Based on files for Hanson & Lee (2019)
%
% Calculates iterated GMM estimator for linear model  
%	y = X*b + e with instrument matrix Z
% Calculates covariance matrix and standard errors robust 
%	to misspecification and clusters
%
% Additionally:
% Calculates the `classic' cluster robust SEs
% Calculates the Windmeijer (2000,2005) corrected SE
%
% Inputs:
%	y       nx1 vector of observations
%	x       nxk matrix of regressors
%	z       nxl matrix of instruments, l>=k (includes exogenous components of x)
%       mem     nx1 vector of cluster membership set mem=(1:n)' for iid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n = size(y,1);
k = size(x,2);
l = size(z,2);
G = length(unique(mem));

% Calculate First Stage
zW      = z'*ww;
zz      = z'*z;
pi1     = (zz)\(zW);

% First Stage Residuals
r       = ww - z*pi1;

% Calculate 2SLS (for comparison)
zx      = z'*x;
zy      = z'*y;
w       = z'*z;
b1      = (zx'/w*zx)\(zx'/w*zy);

% Cluster Index
idx = repmat(mem,1,G)==kron(ones(n,1),unique(mem)');

% Perform the Windmeijer (2019) Weighting with clusters
w = zeros(l,l);
for g = 1:G
    zg  = z(idx(:,g),:);
    eg  = r(idx(:,g));
    zeg = zg'*eg;
    w   = w + zeg*zeg';
end
w = w/n;

% Estimates
b = (zx'/w*zx)\(zx'/w*zy);

% Standard Errors
e = y - x*b;
ze = z.*(e*ones(1,l));
mu = mean(ze)';

if l>k
  J = (mu'/w*mu)*n;
  pv = chi2cdf(J,l-k,'upper');
else
  J = 0;
  pv = 1;
end

if n == G
    ezwze = e.*(z/w*z'*e);
    H = (1/n^2)*(x'*z)/w*(z'*x)-(2/n^3)*(x'*z)/w*(z'*(x.*repmat(ezwze,1,k)));
    Psi = -(1/n)*(z.*repmat(e,1,l))/w*(z'*x)-(1/n)*repmat(((e'*z)/w*z')',1,k).*x+(1/n^2)*repmat(((e'*z)/w*z')',1,k).*((z.*repmat(e.^2,1,l))/w*z'*x);
else
    Hpart = zeros(l,k);
    Psi = zeros(G,k);
    for g = 1:G
        zg = z(idx(:,g),:);
        eg = e(idx(:,g));
        xg = x(idx(:,g),:);
        Hpart = Hpart + (zg'*eg)*(e'*z)/w*(zg'*xg) + (zg'*xg)*((e'*z)/w*(zg'*eg));

        Psi(g,:) = (-(1/n)*(x'*z)/w*(zg'*eg)-(1/n)*(xg'*zg)/w*(z'*e)+(1/n^2)*(x'*z)/w*(zg'*eg)*(eg'*zg)/w*(z'*e))';
    end
    H = (1/n^2)*(x'*z)/w*(z'*x)-(1/n^3)*(x'*z)/w*Hpart;

end

Om = (Psi'*Psi)/n;

* Hanson & Lee misspecification-and-cluster-and-heteroskedasticity robust SE
* note, this is a 2step inefficient GMM, so not exactly sure if correct
V = H\Om/H';
s = sqrt(diag(V/n));
t = b1./s;

* `classic' cluster-and-heteroskedasticity robust SEs
* note, this is a 2step inefficient GMM, so not exactly sure if correct
Q = -zx/n;
V0 = inv(Q'/w*Q);
s0 = sqrt(diag(V0/n));
t0 = b1./s0;

* Windmeijer-corrected cluster-and-heteroskedasticity robust SEs
* note, this is a 2step inefficient GMM, so not exactly sure if correct
Vw = H\(Q'/w*Q)/H';
sw = sqrt(diag(Vw/n));
tw = b1./sw;


