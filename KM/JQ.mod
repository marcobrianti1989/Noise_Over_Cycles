
%% CHECK timing of r

%% Labelling block
% declare endogenous variables
var w c n r b d k mu phip m phi xi z;  

% declare exogenous variables
varexo ez exi; 

% declare parameters
parameters ALPHA TAU THETA DELTA BETA DBAR XIBAR KAPPA RHOZ RHOXI VOLZ VOLXI;


%% Calibration
ALPHA = 1.8834;
TAU = 0.35;     
THETA = 1 - 0.36; 
DELTA = 0.025;
BETA = 0.9825;
XIBAR = log(0.1634);
DBAR = 0;  %must be d_steady state, I am sure there is a way to set it endogenously
KAPPA = 0.146;
RHOZ = .9;
RHOXI = .9;
VOLZ = 1;
VOLXI = 1;


%% Model block
model;
w/c = ALPHA/(1-n); % (1) intratemporal labor consumption

c(+1)/c = BETA*(r-TAU)/(1-TAU); % (2) Euler

w*n + b(-1) - b/r + d - c = 0; % (3) BC 

THETA*exp(z)*n^(THETA -1)*k(-1)^(1-THETA) = w*(1/(1-mu*phip)); % (4) labor demand

m*(1 - DELTA + (1 - mu(+1)*phip(+1))*(1-THETA)*exp(z)*n^THETA*k(-1)^(-THETA)) + exp(xi)*mu*phip = 1; % (5) firm foc for kp
 
r*m + exp(xi)*mu+phip*(r*(1-TAU)/(r-TAU)) = 1; % (6) 

(1-DELTA)*k(-1) + exp(z)*n^THETA*k(-1)^(1-THETA)-w*n-b(-1)+b/r-k-phi = 0; %(7)

exp(xi)*(k-b*(1-TAU)/(r-TAU)) = exp(z)*n^THETA*k(-1)^(1-THETA); %(8)

m = BETA*(c/c(+1))*(phip/phip(+1)); 

phi = d + KAPPA*(d - DBAR)^2; 

phip = 1 + 2*KAPPA*(d-DBAR);

exp(z) = exp(RHOZ*z(-1) + ez); 

exp(xi) = exp(RHOXI*xi(-1) + XIBAR*(1-RHOXI) + exi);

end;

%% Initialization block
%initval;

%end;

steady;
check;
%% Random shocks block
shocks;
var ez; stderr VOLZ ;
var exi; stderr VOLXI ;
var ez,exi = 0;
end;


%% Solution and property block
steady;
check;

stoch_simul(periods=200,order=1,irf=25) k; % order is the Taylor approximation; %irf is the # periods 

%options_.noprint = 1;

