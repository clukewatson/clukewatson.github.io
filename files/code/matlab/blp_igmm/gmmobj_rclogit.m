function [f, g] = gmmobj_rclogit( THETANL )
%{
THETANL = thetaNL_init;
MU 	= mu_init;
%}

global invA thetaL X Z fcnevals gmmresid MU PSI;

disp("coefs:[alpha, vec(sigma)]")
disp([THETANL])

% Find the Mean Value That Rationalized Trial Parameters
MU     = meanvalueCM_mix(THETANL);

% Indicators for MU Output
% 1 = if there is a NaN value
% 2 = if there is an infinite value
% 3 = if there is either NaN or Infinite Value
cond1 = double( max( isnan( MU ) ) ) ;
cond2 = double( min( isfinite( MU ) ) == 0 );
cond3 = cond1 + cond2;

% *********************************************************************
% Calculate GMM ObjValue & Gradient
% *********************************************************************

if cond3 >=1

    f=NaN;
    g=NaN;
    fcnevals=fcnevals+1;
    disp("Mean Val NA")

else

% Linear Parameters
    A = X'*(Z*(invA*(Z'*X)));
    B = X'*(Z*(invA*(Z'*MU)));
    thetaL   = A \ B ;

% Use residual to form non-linear parameter moments
    gmmresid = MU - (X*thetaL);
    temp3    = (gmmresid'*Z);

% GMM Objective Function Value
    f        = temp3*(invA*temp3');
    disp("f done")
    disp([f])

% GMM Objective Function Gradient
    g        = gradient_rclogit( MU , PSI );

% Update Evaluation Counter
    fcnevals = fcnevals+1;

% Export Trial Non-linear Parameters
    dlmwrite('thetaNL_out.dat',THETANL)

end
