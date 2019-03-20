% Brianti, Cormun January 2019 - Buaudry, Galizia, Portier (2019, AER)

%%%%%%% Defining Variables %%%%%%%

var 
                              
f                      % log of financial fragility
b                      % intraperiod bond
c                      % Consumption
k                      % capital
n                      % hours workerd
w                      % wage  
y                      % output
i                      % investment
r                      % interest rate
logz                   % stochastic productivity
logF;                  % stochastic fixed cost

%%%%% Aggregate Productivity Shock %%%%%

varexo epsz epsF;  

% Parameterization
parameters
Z EPS ETA FF ALP BET GAM DEL SIG PSI CHI RHOZ RHOF RHOf; %XX

Z          = 1;      % mean of productivity of all the firms 
ETA        = 0.1;   % elasticity of financial fragility to past output
EPS        = 0.8;    % percentage of capital recover if default
FF         = 1;      % fixed cost of recovering capital if default
ALP        = 2/3;    % Convexity of Production Function
BET        = 0.99;   % Discount Factor
GAM        = 0.8;    % habit consumption formation
DEL        = 0.05;   % Capital depreciation 
SIG        = 2;      % IES
PSI        = 1;      % multiplicative disutility of labor
CHI        = 4;      % Frish elasticity
RHOZ       = 0.95;   % persistence of productivity shock
RHOF       = 0.2;   % persistence of sentiment shock 
RHOf       = 0.5;    % persistence of financial fragility

%XX        = 1;      % fixed cost of production

% Defining functional forms and derivatives

model; 

b = (Z*exp(logz)*k^ALP*n^(1-ALP))*EPS - FF; %*exp(f)

f = RHOf*f(-1) - ETA*log(k(-1));

b = w*n + r*k;

PSI*n^CHI = w;

c + i = y;

ALP*n*w = (1-ALP)*r*k;

1 = BET*((c - GAM*c(-1))/(c(+1)- GAM*c))^SIG*(1+r(+1)-DEL)*exp(logF);

y = Z*exp(logz)*k^ALP*n^(1-ALP);

k(+1) = (1-DEL)*k + i;

logz       = RHOZ*logz(-1) + epsz;  

logF       = RHOF*logF(-1) - epsF;            



end;

rss    = 1/BET - 1 + DEL;
css    = 0.0172;
nss    = 0.3704;
kss    = 1/rss*ALP/(1-ALP)*PSI*nss^(1+CHI);
wss    = PSI*nss^CHI;
yss    = Z*kss^ALP*nss^(1-ALP);
iss    = kss*DEL;
bss    = wss*nss + rss*kss;
logzss = 0;
logFss = 0;
fss    = 0;
                         
%%%%% Initialization %%%%%
initval;

k                = kss; 
c                = css;
r                = rss;
y                = yss;
i                = iss; 
b                = bss;
n                = nss;
w                = wss; 
logz             = 0;
logF             = 0;
f                = 0;

end;


shocks;
  var epsz     = 1;
var epsF = 1;
end;

steady;
check;

stoch_simul(irf=60, order=1) y c i k r w n b;


%stoch_simul(periods=100000, hp_filter = 1600, order=2,ar = 0, nofunctions,nograph,nodecomposition,nocorr) jobphim logz logR logU logV logFPHIC;



