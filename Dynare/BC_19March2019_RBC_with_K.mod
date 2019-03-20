% Brianti, Cormun January 2019 - Buaudry, Galizia, Portier (2019, AER)

%%%%%%% Defining Variables %%%%%%%

var 
                              
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
Z ALP BET GAM DEL SIG CHI RHOZ RHOF;

Z          = 2;      % mean of productivity of all the firms 
ALP        = 2/3;    % Convexity of Production Function
BET        = 0.99;   % Discount Factor
GAM        = 0.7;    % habit consumption formation
DEL        = 0.05;   % Capital depreciation 
SIG        = 2;      % IES
RHOZ       = 0.7;    % persistence of productivity shock
RHOF       = 0.7;    % persistence of sentiment shock 
CHI        = 15;      % procyclical marginal cost

% Defining functional forms and derivatives

model; 

r = exp(logz)*ALP*k^(ALP-1);

c + i = y;

1 = BET*((c - GAM*c(-1))/(c(+1)- GAM*c))^SIG*(1+r(+1)-DEL)*exp(logF);

y = 1/2*Z*exp(logz)*k^ALP;

k = (1-DEL)*k(-1) + i(-1);

logz = RHOZ*logz(-1) + epsz;  

logF = RHOF*logF(-1) - epsF;            

end;

rss    = 1/BET - 1 + DEL;
kss    = 133.1092;
css    = 18.9591;
yss    = 1/2*Z*kss^ALP;
iss    = kss*DEL;
                         
%%%%% Initialization %%%%%
initval;

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

stoch_simul(irf=20, order=1) k r c y i;


%stoch_simul(periods=100000, hp_filter = 1600, order=2,ar = 0, nofunctions,nograph,nodecomposition,nocorr) jobphim logz logR logU logV logFPHIC;



