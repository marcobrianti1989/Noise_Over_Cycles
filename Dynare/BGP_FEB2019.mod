% Brianti, Cormun January 2019 - Buaudry, Galizia, Portier (2019, AER)

%%%%%%% Defining Variables %%%%%%%

var 
c
rp                     % (1) Risk Premium
e                      % (2) Employment Rate 
x                      % (3) Durable Goods Stock 
logthet                % (4) Log Technology
logzeta;               % (5) Log Preferences

%%%%% Aggregate Productivity Shock %%%%%

varexo eps_thet eps_zeta;  

% Parameterization
parameters
OMEG GAM PSI PHIE PHI PHIBIG S ALP DEL THET BET RHO_THET SIGMA_THET RHO_ZETA SIGMA_ZETA; 

OMEG       = 0.2408; % CRRA Parameter
GAM        = 0.5;     %0.5876; % Habit
PSI        = 0.2994; % 1 - depreciation of durables in the first period 
PHIE       = 0.0467; % Taylor Rule Parameter
PHI        = 0.95; %0.8827; % Probability of recovering impaired loans
PHIBIG     = 0.0458; % Recovery Cost
S          = 1;      % Other Tech Parameter
ALP        = 2/3;    % Convexity of Production Function
DEL        = 0.05;   % Durable Goods Depreciation Rate 
BET        = 0.99;   % It could be WRONG! Discount Factor
<<<<<<< HEAD
THET       = 1/(BET*(1 - 0.0583)^PHIE)*.9976; % It should be set to have a ss unemployment rate of 0.0583
RHO_THET   = 0.85;    % Persistence of tech shock
SIGMA_THET = 1;      % SD of tech shock
RHO_ZETA   = 0.74;    % Persistence of Preference Shocks
SIGMA_ZETA = .01;      % SD of preference shock  
=======
THET       = 1/(BET*(1 - 0.0583)^PHIE); % It should be set to have a ss unemployment rate of 0.0583
RHO_THET   = 0.9;    % Persistence of tech shock
SIGMA_THET = 1;      % SD of tech shock
RHO_ZETA   = 0.4;    % Persistence of Preference Shocks
SIGMA_ZETA = 1;      % SD of preference shock  
>>>>>>> e006909e2d0d57d43d122b1a743af23eb2277219

% Defining functionals forms and derivatives

model; 

rp = (1 + (1 - e)*PHI*PHIBIG )/(e + (1 - e)*PHI);   %(1)

x = (1 - DEL)*x(-1) + PSI*exp(logthet(-1))*e(-1)^ALP;          %(2)

logthet       = RHO_THET*logthet(-1) + SIGMA_THET*eps_thet;  %(3) 

logzeta       = RHO_ZETA*logzeta(-1) - SIGMA_ZETA*eps_zeta;  %(4) 

(S*(x + exp(logthet)*e^ALP) - GAM*S*(x(-1) + exp(logthet(-1))*e(-1)^ALP))^(-OMEG) = BET*THET*exp(logzeta)/exp(0*logzeta(-1))*(e + (1-e)*PHI)*rp*(S*(x(+1) + exp(logthet(+1))*e(+1)^ALP) - GAM*S*(x + exp(logthet)*e^ALP))^(-OMEG)*e(+1)^PHIE; %(5)

c = x + exp(logthet)*e^ALP; %(5) consumption
end;

logthetss        = 0;
logzetass        = 0;
ess              = 0.777146777146777;
rpss             = 1.03609362949245;
xss              = 5.06155342500016;
css              = xss + ess^ALP;                                
%%%%% Initialization %%%%%
initval;

logthet              = logthetss;  
logzeta              = logzetass;
e                    = ess;
rp                   = rpss;
x                    = xss;
c                    = css;
end;


shocks;
  var eps_thet     = 1;
  var eps_zeta    = 1;
  var eps_thet, eps_zeta = 0;
end;

steady;
check;

<<<<<<< HEAD
stoch_simul(irf=30, order=1) logthet logzeta e x c;
=======
stoch_simul(irf=20, order=1) logthet logzeta e x rp;
>>>>>>> e006909e2d0d57d43d122b1a743af23eb2277219

%stoch_simul(periods=100000, hp_filter = 1600, order=2,ar = 0, nofunctions,nograph,nodecomposition,nocorr) jobphim logz logR logU logV logFPHIC;



