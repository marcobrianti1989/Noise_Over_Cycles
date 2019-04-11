% Kiyotaki Moore model with exogenous changes in supply of land and sentiment shocks 
% Model delivers predictions on real down-payments (could also be interpreted as easiness to access credit)
% derive the reduced form to understand whether is consistent with our initial intuition 
% derive implications for endogenous tfp 
% need to show that the parametrization is such that the constraint always binds
% add sentiment shock check azariadis et al
% what's the externality?
% understand the importance of the lack of uncertainty in this model
% consumption and investment tradeoff in response to demand shocks 
% downpayment in the data
% review capital accumulation and understand whether they already relax the assumption of fixed land. Write a summary on the model with capital 
% See Krishnamurthy  JET. Iacoviello 2005 AER: KM with sticky prices, codes online
% Ireland 2004 on tech shocks Restat


%% Labelling block
% declare endogenous variables
var k q b z kbar s dp yg y c i cg cons;  

% declare exogenous variables
varexo ez es; 

% declare parameters
parameters PI LAMBDA PHI R A M ALPHA VOLZ RHOZ RHOS VOLS C;


%% Calibration
A   = .5; %1
PHI = 1.4; %1.3; %20
R   = 1.025; %1.01
LAMBDA = 0.9875;
PI = 0.60 ;  
ALPHA  = .5/3;
M = 1;
VOLZ = .1;
VOLS = 0.01;
RHOZ = .92;
RHOS = .92;
C = .5;

%% Model block
model;

k =  (PI/(q + PHI - q(+1)/R))*( (A*exp(z)+LAMBDA*PHI+q)*k(-1) - R*b) + (1-PI)*LAMBDA*k(-1);

b*exp(s)    =  R*b(-1) + q*(k - k(-1)) + PHI*(k-LAMBDA*k(-1)) - A*exp(z)*k(-1);

ALPHA*exp(z(+1))*(((kbar - k)/M)^(ALPHA-1))/R   =  q - q(+1)/R; 

kbar =2; %*exp(10*z);  

yg = M^(1-ALPHA)*exp(z)*(kbar(-1)-k(-1))^ALPHA; %gatherers'aggregate output

y = yg + (A+C)*exp(z)*k(-1); %aggregate output

c = yg + A*exp(z)*k(-1)-PHI*(k-LAMBDA*k(-1)) + C*exp(z)*k(-1); %aggregate consumption

cg = c -  C*exp(z)*k(-1);

exp(z) = exp(RHOZ*z(-1) + ez);

exp(s) = exp(RHOS*s(-1) + es);

dp = q - q(+1)/R;

i = PHI*(k-LAMBDA*k(-1));

cons = q*k(-1) - R*b(-1);

end;

%% Initialization block
initval;
q = (R/(R-1))*(PI*A - (1-LAMBDA)*(1- R + PI*R)*PHI )/(LAMBDA*PI + (1-LAMBDA)*(1- R + PI*R));
b  =(A + LAMBDA*PHI - PHI)/(R-1) * k;
z  = 0;
kbar = 2;
k  = kbar - M*(R*q/ALPHA - q/ALPHA)^(1/(ALPHA-1));
s = 0;
yg = M^(1-ALPHA)*(kbar-k)^ALPHA; 
y = yg + (A+C)*k; 
c = yg + A*k-PHI*(k-LAMBDA*k) + C*k; 
dp = q*(1-1/R);
i = PHI*(k-LAMBDA*k);
cg = c -  C*k;
cons = q*k - R*b;
end;

steady;
check;
%% Random shocks block
shocks;
var ez; stderr VOLZ ;
var es; stderr VOLS ;
var ez,es = 0;
end;


%% Solution and property block
steady;
check;

stoch_simul(periods=200,order=1,irf=25) k b q z s i  y c dp cons; % order is the Taylor approximation; %irf is the # periods 

%options_.noprint = 1;

