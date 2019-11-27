function [PSI_ji] = psifunc(THETANL)
%{
THETANL = thetaNL_init;
thNL_p = -9.0;
%}

global vi yi Xc P kc ;

thNL_p 	= THETANL(1);
thNL_v 	= THETANL(2:end);

vv 	= Xc(:,1:kc)*(vi' .* thNL_v');
pp 	= (thNL_p)*(P./yi');
PSI_ji 	= vv + pp;

end


