
clear all

%**************************************************************************
%**************************************************************************
% PART ONE
%**************************************************************************
%**************************************************************************


% *********************************************************************
% Define globals
% *********************************************************************
global ns k kc K J X Xc Z P SH_j ...
	vi yi PSI MU mu_init ...
	invA thetaL gmmresid ...     	
	step_max step_min ...
     	mymaxfunevals fcnevals ;

% *********************************************************************
% Program Parameters
% *********************************************************************

% Number of Simulated Consumers
ns  = 3000;

% Optimization Parameters
mytolx        = 1e-6;
mytolfun      = 1e-6;
mymaxiters    = 5*10^5;
mymaxfunevals = 5000;
step_max      = 1;
step_min      = 1;

    
% *********************************************************************
% Load data
% *********************************************************************

% File Path
cd '.../rclogit/';
EST_DATA = readtable("est_data.csv");


% *********************************************************************
% Observations, Shares, Initial Mean Value Guess, Price
% *********************************************************************
%
% Number of Buildings
prod_id = EST_DATA(:,'id');
J = size(prod_id,1);

% Variables
% unit price, price-trimmed
P = (table2array(EST_DATA(:,'Price')));

% True Share
SH_j = table2array(T(:,'sh_HH_j'));

% Initial Mean Value Guess = Linear Utility Difference
psi_j = table2array(T(:,'psi_HH_j'));

% *********************************************************************
% Exogenous Variables in Utility Function 
% *********************************************************************
%
% exogenous variables : includes constant
Xexog_continuous = table2array(T(:,{'X0' 'X1' 'X2' 'X3' 'X4' 'X5'}));
Xexog_catagorical = table2array(T(:,{'D1' 'D2'}));
sdXsq_c = std(Xexog_continuous(:,2:end));
Xexog_continuous(:,2:end) = Xexog_continuous(:,2:end) ./ sdXsq_c;

% Check distribution of variables
% [mean(Xexog_continuous), std(Xexog_continuous), min(Xexog_continuous), max(Xexog_continuous)]
% [mean(Xexog_catagorical), std(Xexog_catagorical)]

% Fixed Effect Variable
FE = table2array(T(:,'FEVAR'));
[~, ~, FE_idx] = unique(TR, 'rows');
iFE = dummyvar(FE_idx);
kFE = size(iFE,2);

% check is they make sense
FEsum = sum(iFE)';
FEsum = FEsum(FE_idx)';

% Exogenous Variables Vector
X = [Xexog_continuous Xexog_catagorical iFE(:,2: kFE)];

% Only Continuous X get random coefs
% Xc = table2array(T(:,{'X0' 'X1' 'X2' 'X3' 'X5'}));
% sdXc = std(Xc(:,2:end));
% Xc(:,2:end) = Xc(:,2:end) ./ sdXc;
Xc = Xexog_continuous ; 


% *********************************************************************
% Exogenous Variables for Price Equation
% *********************************************************************
%
%%% Identification Conditions
% number of continuous X
kc   = size(Xc,2)
% Total number of X
k   = size(X,2)

% Total Number of Model Parameters
% price, nest, mean_X, sd_X
K   = 1 + 1 + k + kc*(R + 1)
% number of coefs the instruments need to identify
% alpha, rho, sigma0-sigma7
Kz  = 1 + 1 + kc*(R + 1)

% Cost Shifting Instruments
Z_cost_shift = table2array(T(:,{'Z1' 'Z2'}));
                     
% Substitution Shifting Instruments
% -- made elsewhere
% -- usually at least one for each continuous variable
% -- here there is also a predicted price and IVs for Indicators
Z_sub_shift = table2array(T(:,{ ...
			'iv_X0' 'iv_X1' 'iv_X2' 'iv_X3' 'iv_X4' ... 
			'iv_X5' 'iv_X6' 'iv_D1' 'iv_D2' 'iv_Dq_P1' }));

% Scale by Standard Deviation like other continuous variables
tZ =   [Z_cost_shift Z_sub_shift];
sdZsq = std(tZ);
tZ = tZ ./ sdZsq;

% Make Instrument Vector
Z = [X tZ];


%**************************************************************************
%**************************************************************************
% PART TWO
%**************************************************************************
%**************************************************************************


% *********************************************************************
% Use Halton Draws 
% *********************************************************************

% Halton Sequence for Random Draws
skip_length = 10000;
leap_length = 52;
p = haltonset(kc+1,'Skip',skip_length,'Leap',leap_length);
p = scramble(p,'RR2');

% Creates (ns)x(K+3) matrix of (0,1) points
u_halton = p(1:ns,:);
clear p;

% Check that they look like uniform draws
%mean(u_halton)
%var(u_halton)
%boxplot(u_halton)


% *********************************************************************
% Draws for Income
% *********************************************************************
%
% y ~ lognormal( mean(l og[y] ) , var( log[y] ) )
% --> log[y]  ~ normal( mean( log[y] ) , var( log[y] ) )
% Let
% mean( log[y] ) == MU
% var( log[y] ) == VAR
%
% Note, if income variable is truncated by bottom and top coding
% 	then for the summary statistics:
%	mu_y is biased up and var_y is ambiguously biased.
% With assumptions, one can model this, 
% 	but no one in the literature seems to care about this

z_yi 	= norminv(u_halton(:,1));
yi 	= exp( MU  + ( VAR.*z_yi ) );
clear z_yi;

% Check that this looks correct
% boxplot(yi)
% mean(yi)
% median(yi)

% *********************************************************************
% Draws for the Sigma_v
% *********************************************************************
% 
vi = norminv(u_halton(:,2:end)); % Nxk

% Check that this looks correct
% mean(vi)
% var(vi)
% boxplot(vi)

clear u_halton;

%**************************************************************************
%**************************************************************************
% PART THREE
%**************************************************************************
%**************************************************************************


% *********************************************************************
% Optimization Parameters and Initial Values
% *********************************************************************

% Initial THETANL
thetaNL_init = [ -10.0 ones(1,kc) ];

% initialize counter of function evaluations
fcnevals    = 0;

% Inital Weight Matrix
ZZ      = Z'*Z;
invA    = ZZ\eye(size(ZZ));

% initialize MU
mu_init     = psi_j;
MU          = mu_init;
PSI         = [];

% fmincon constraints
A   = [];
b   = [];
Aeq = [];
beq = [];
lb  = [  -100.0,  0*ones(1,kc) ];
ub  = [  -  1.0, 10*ones(1,kc) ];
nonlcon = [];
% fmincon options
options = optimoptions( 'fmincon' , ...
		 	'SpecifyObjectiveGradient' , true, ...
                        'MaxFunEvals' , mymaxfunevals, ...
                        'TolFun'  , mytolfun, ...
                        'TolX'    , mytolx, ...
                        'MaxIter' , mymaxiters, ...
                        'Display' , 'iter' ) ;

% Containers
fval_iters      = zeros(1,100);
mean_ope_iter   = zeros(1,100);
thetaNL_iters   = zeros(kc+1,100);
normTH_iter     = zeros(1,100);
normWW_iter     = zeros(1,100);

% begin
normTH      = 1;
i           = 1;
thetaNL0    = thetaNL_init;
invA0       = invA;

% *********************************************************************
% LOOP
% *********************************************************************

while (normTH > 0.01 && i<=100)

% ****************
% GMM Optim fmincon
% ****************
[ thetaNL1 , fval , exit_info , tmp ] = ...
                    fmincon( 'gmmobj_rclogit' , ...
                              thetaNL_init, ...
                              A,b,Aeq,beq,lb,ub,nonlcon,options );
                            
% ****************
% Output
% ****************
fval
exit_info
tmp

disp(thetaNL1')

thetaNL_iters(:,i)  = thetaNL1';
fval_iters(i)       = fval;

PSI         = psifunc( thetaNL1 );
MU          = meanvalueCM( thetaNL1 );
[SJ]        = ind_share( MU , PSI );
mShj        = mean( SJ , 2 );
alpha       = thetaNL1( 1 );
ope_j1      = (alpha)*(P./mShj).*(mean((1./yi)'.*SJ.*(1 - SJ),2));
clear SJ ;

disp(["Own Price Elasticity Parameters"; ...
    "min max mean sd median"])
fprintf('%.7f\n',[min(ope_j1) max(ope_j1)  ...
    mean(ope_j1) sqrt(var(ope_j1)) median(ope_j1)])

mean_ope_iter(i) = mean(ope_j1);

% ****************
% Reset Weight Matrix 
% for next iteration
% ****************
Am      = X'*(Z*(invA*(Z'*X)));
Bm      = X'*(Z*(invA*(Z'*MU)));
thetaL  = Am \ Bm ;
gmmresid = MU - X*thetaL;
ZuuZ    = (Z' * (diag(gmmresid.^2) * Z));
tmp     = (ZuuZ + ZuuZ')/2;
invA1   = tmp\eye(size(tmp));

% ****************
% Evaluate iteration
% ****************

% Norms
normTH              = norm((thetaNL1 - thetaNL0));
normWW              = norm((invA1 - invA0));
normTH_iter(i)      = normTH;
normWW_iter(i)      = normWW;

disp("GM:normTH, normWW, mean ope, alpha, fval")
%disp([i, normTH, normWW])
disp([normTH_iter(1:i)' normWW_iter(1:i)' mean_ope_iter(1:i)' thetaNL_iters(1,1:i)' fval_iters(1:i)'])

disp([thetaNL_iters(:,1:i)])

% Set the "OLD"    
thetaNL0    = thetaNL1;
invA0       = invA1;

% Reset Weight Matrix
invA        = invA1;

% counter
i = i+1;

% ****************
% Export Iteration Values
% ****************

csvwrite("fval_iters.csv",fval_iters) 
csvwrite("mean_ope_iter.csv",mean_ope_iter) 
csvwrite("thetaNL_iters.csv",thetaNL_iters) 
csvwrite("normTH_iter.csv",normTH_iter) 
csvwrite("normWW_iter.csv",normWW_iter) 

end

% *********************************************************************
% Loop complete
% *********************************************************************


%**************************************************************************
%**************************************************************************
% PART FOUR
%**************************************************************************
%**************************************************************************


thetaNLopt = thetaNL1;

% *********************************************************************
% Final Elasticity
% *********************************************************************
%
PSI         = psifunc(thetaNLopt);
MU          = meanvalueCM(thetaNLopt);
[SJ]        = ind_share(MU,PSI);
mShj        = mean(SJ,2);
alpha       = thetaNLopt(1);
ope_j1      = (alpha)*(P./mShj).*(mean((1./yi)'.*SJ.*(1 - SJ),2));
clear SJ ;

disp("If it finishes")
disp("Estimates of NL Parameters")
disp([thetaNLopt' thetaNL_iters(:,2) thetaNL_iters(:,1) thetaNL_init'])
disp(["Own Price Elasticity Parameters"; ...
    "min max mean sd median"])
fprintf('%.7f\n',[min(ope_j1) max(ope_j1)  ...
    mean(ope_j1) sqrt(var(ope_j1)) median(ope_j1)])

% *********************************************************************
% Standard Errors
% *********************************************************************
%
Am      = X'*(Z*(invA*(Z'*X)));
Bm      = X'*(Z*(invA*(Z'*MU)));
thetaL  = Am \ Bm ;
gmmresid = MU - X*thetaL;
dj      = jacob_rclogit(MU,PSI);
ZX      = ([(-1*X) dj]'*Z);
ZuuZ    = (Z' * diag(gmmresid.^2) * Z);

BRD     = ZX*(invA*ZX'); 
iBRD    = BRD\eye(size(BRD));
Meat    = ZX*(invA*(ZuuZ*(invA*ZX')));
V       = iBRD*(Meat*iBRD);

se      = (full(sqrt(diag(V))));
theta   = [thetaL; thetaNLopt'];
t       = theta ./ se; 

% *********************************************************************
% Export Results
% *********************************************************************
%
gmm_out = [theta, se, t];
dlmwrite('gmm_out_opt.csv',gmm_out,'delimiter',',') 
csvwrite("gmm_out2.csv",gmm_out) 

ope_out = [bbl , array2table(P), array2table(ope_j1)];
writetable(ope_out,'ope_out_opt.csv','delimiter',',') 

% 
