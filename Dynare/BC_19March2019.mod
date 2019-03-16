% Brianti, Cormun January 2019 - Buaudry, Galizia, Portier (2019, AER)

%%%%%%% Defining Variables %%%%%%%

var 
                              
rho                    % Probability of default
y                      % Production
logz                   % stochastic productivity
logF                   % stochastic fixed cost
n;                     % labor

%%%%% Aggregate Productivity Shock %%%%%

varexo epsz epsF;  

% Parameterization
parameters
Z X EPS F ALP BET PSI CHI RHOZ RHOF;

Z          = 2;      % mean of productivity of all the firms 
X          = 0.01;   % fixed cost of production
EPS        = 0.7;    % percentage of capital recover if default
F          = 1;      % fixed cost of recovering capital if default
ALP        = 2/3;    % Convexity of Production Function
BET        = 0.99;   % Discount Factor
PSI        = 0.1;    % linear consumption disutility
CHI        = 4;      % Frish
RHOZ       = 0;   % persistence of productivity shock
RHOF       = 0;   % persistence of sentiment shock 

% Defining functional forms and derivatives

model; 

PSI*(1 - rho^2)*n^(1+CHI) = EPS*Z*exp(logz)*rho*n^ALP - F*exp(logF);

rho*n^ALP - PSI*(1 - rho^2)*n^(1+CHI) - X/(1/2*Z*(1-rho(-1)^2)*n(-1)^ALP) = 0;

y = (1 - rho)*n^ALP;

logz       = RHOZ*logz(-1) + epsz;  

logF       = RHOF*logF(-1) - epsF;            



end;

rhoss = 0.3804  
nss   = 0.7396
                         
%%%%% Initialization %%%%%
initval;

rho              = rhoss;  
n                = nss; 

end;


shocks;
  var epsz     = 1;
var epsF = 1;
end;

steady;
check;

stoch_simul(irf=20, order=1) n rho y;


%stoch_simul(periods=100000, hp_filter = 1600, order=2,ar = 0, nofunctions,nograph,nodecomposition,nocorr) jobphim logz logR logU logV logFPHIC;



