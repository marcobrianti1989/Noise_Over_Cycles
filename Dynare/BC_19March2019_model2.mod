% Brianti, Cormun January 2019 - Buaudry, Galizia, Portier (2019, AER)

%%%%%%% Defining Variables %%%%%%%

var 
                              
rho                    % Probability of default
c                      % Consumption
k                      % capital
y                      % output
i                      % investment
r                      % interest rate
logz                   % stochastic productivity
logF;                  % stochastic fixed cost

%%%%% Aggregate Productivity Shock %%%%%

varexo epsz epsF;  

% Parameterization
parameters
Z XX EPS FF ALP BET GAM DEL SIG RHOZ RHOF;

Z          = 2;      % mean of productivity of all the firms 
XX         = 1;   % fixed cost of production
EPS        = 0.8;    % percentage of capital recover if default
FF         = 0.2;      % fixed cost of recovering capital if default
ALP        = 2/3;    % Convexity of Production Function
BET        = 0.99;   % Discount Factor
GAM        = 0.7;    % habit consumption formation
DEL        = 0.05;   % Capital depreciation 
SIG        = 2;      % IES
RHOZ       = 0.5;      % persistence of productivity shock
RHOF       = 0.5;      % persistence of sentiment shock 

% Defining functional forms and derivatives

model; 

r*k(-1) = EPS*Z*exp(logz)*rho*k(-1)^ALP - FF;

exp(logz)*Z*rho*k(-1)^ALP - r*k(-1) = XX;

c + i = y  - (1 - rho)*XX;

1 = BET*((c - GAM*c(-1))/(c(+1)- GAM*c))^SIG*(1+r(+1)-DEL)*exp(logF);

y = 1/2*Z*exp(logz)*(1 - rho^2)*k^ALP;

k = (1-DEL)*k(-1) + i;

logz       = RHOZ*logz(-1) + epsz;  

logF       = RHOF*logF(-1) - epsF;            



end;

rhoss  = 0.1726;
kss    =  133.1092;
css    = 18.9591;
rss    = 1/BET - 1 + DEL;
yss    = 1/2*Z*(1 - rhoss^2)*kss^ALP;
iss    = kss*DEL;
                         
%%%%% Initialization %%%%%
initval;

rho              = rhoss;  
k                = kss; 
c                = css;
r                = rss;
y                = yss;
i                = iss; 

end;


shocks;
  var epsz     = 1;
var epsF = 1;
end;

steady;
check;

stoch_simul(irf=60, order=1) rho k r c y i;


%stoch_simul(periods=100000, hp_filter = 1600, order=2,ar = 0, nofunctions,nograph,nodecomposition,nocorr) jobphim logz logR logU logV logFPHIC;



