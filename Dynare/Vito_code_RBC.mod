% RBC Model With Isoelastic Utility (PHI = 0.25)
 
var c g h i k r w y z;
varexo ez eg ;
 
parameters ALPHA, BETA, DELTA, PHI, RHOG, RHOZ, SG, SC, VOLG, VOLZ;
 
ALPHA = 1/3;
BETA  = 0.99;
DELTA = 0.025;
PHI   = 0.25;
RHOG  = 0.90;
RHOZ  = 0.99;
SG    = 0.20;
SC    = 0.65;
VOLG  = SG^(-1)/100;
VOLZ  = (1-ALPHA)^(-1)/100;
 
 

model(linear);
k = (1-DELTA)*k(-1) + DELTA*i;
y = SC*c + SG*g + (1-SC-SG)*i;
y = ALPHA*k(-1) + (1-ALPHA)*z + (1-ALPHA)*h;
w = ALPHA*k(-1) + (1-ALPHA)*z - ALPHA*h;
r = -(1-ALPHA)*k(-1) + (1-ALPHA)*z + (1-ALPHA)*h;
c + PHI* h = w;
c(+1) - c = (1-BETA+BETA*DELTA)*r(+1);
z = RHOZ*z(-1) + ez;
g = RHOG*g(-1) + eg;
end;
 
initval;
c = 0;
g = 0;
h = 0;
i = 0;
k = 0;
r = 0;
w = 0;
y = 0;
z = 0;
end;
 
shocks;
var ez; stderr VOLZ ;
var eg; stderr VOLG ;
var ez,eg = 0;
end;
 
steady;
check;
 
%options_.noprint = 1;
 
stoch_simul(periods=200,order=1,irf=40,nocorr) y c i h w k r z g;