function f = meanvalueCM_mix(THETANL)
%{
This function performs the contraction mapping and returns mean value mu_j.
This code refers to this as mval for mean value.
This code uses the SQUAREM routine.
ALWAYS starts from logit Mu (mu_init)

THETANL = thetaNL_init;
THETANL = thetaNL0;

MUold = mu_init;
step_min =1;
step_max =1; 
%}

global SH_j PSI step_max step_min mu_init;

% *********************************************************************
% Parts for Nested Fixed Point
% *********************************************************************

% Individual Heterogeneity
PSI = psifunc(THETANL);

% Log Market Share (should global this)
lsh = log(SH_j);

% Start Point : Logit MU
MUiter = mu_init;

% *********************************************************************
% Nested Fixed Point
% *********************************************************************

% 1) Try SQUAREM
% 2) If SQUAREM fails, then use BLP Inversion

% Params
step_factor =4;
step_minT = step_min;
step_maxT = step_max; 
norm = 1;
i = 0;
tic

while (norm > 1e-12 && i<=1000)
    
          A  = MUiter + (lsh - log(mrktshare(MUiter,PSI)));
          
          % Next, Next Iteration
          B = A + (lsh - log(mrktshare(A,PSI)));
          
          % Resids
          r1 = (A - MUiter);
          r2 = (B - A) - r1;
          
          % Create Step Length
          a1 = - sqrt((r1'*r1)/(r2'*r2));
          a2 = -1*(max(step_minT, min(step_maxT, -1*a1)));
          
          % Acceleration
          C = MUiter - 2*a2.*r1 + (a2.^2).*r2;
          D = C + (lsh - log(mrktshare(C,PSI)));
          
          normC = max(abs(D-C));
          
          % Evaluate Step
          if abs(a2)==step_maxT 
              step_maxT= step_factor*step_maxT;
          end
          
          % Evaluate Iteration
          norm = normC;
          meanD = mean(D);
          MUiter = D;
          i = i + 1;

         %disp([i, mean(D),step_maxT])
end

tmp1 = isfinite(norm);
tmp2 = step_maxT>64;

if (tmp1==0) || (tmp2==1)
    disp("SQEM FAILED")

    % Tolerances
    norm = 1;
    i = 0;

    % Start Value
    MUiter = mu_init;
    
    while (norm > 1e-12 && i<=2000)

	% Step
	MUnew = MUiter + (lsh - log(mrktshare(MUiter,PSI)));

	% Evaluate  
	t = abs(MUnew-MUiter);
	norm = max(t);

	% Reset
	MUiter = MUnew;
	i = i + 1;

	% Report
	%[i , norm, mean(MUiter)]

    end      

end 
elapsedTime = toc;
disp("CM: i, min, mean, step")
disp([i, elapsedTime / 60, mean(D),step_maxT])

f = (D);

end
