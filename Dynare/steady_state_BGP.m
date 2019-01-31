
% Steady State Parameters
OMEG       = 0.2408; % CRRA Parameter
GAM        = 0.5876; % Habit
PSI        = 0.2994; % One minus initial debt
PHIE       = 0.0467; % Taylor Rule Parameter
PHI        = 0.8827; % Debt Baking
PHIBIG     = 0.0458; % Recovery Cost
S          = 1;      % Other Tech Parameter
ALP        = 2/3;    % Convexity of Production Function
DEL        = 0.05;   % Durable Goods Depreciation Rate 
BET        = 0.99;   % It could be WRONG! Discount Factor
THET       = 1/(BET*(1-0.0583)^PHIE); % It should be set to have a ss unemployment rate of 0.0583
RHO_THET   = 0.9;    % Persistence of tech shock
SIGMA_THET = 1;      % SD of tech shock
RHO_ZETA   = 0.9;    % Persistence of Preference Shocks
SIGMA_ZETA = 1;      % SD of preference shock  


e = linspace(0,1,1000000);
check = (  (1 + (1 - e).*PHI*PHIBIG).*BET.*THET.*e.^PHIE - ones(1,length(e))  ).^2;

[zerocheck, loc_estar] = min(check);

estar = e(loc_estar);

rp = (1 + (1 - estar)*PHI*PHIBIG)/(estar + (1 - estar)*PHI);

x = PSI/DEL*estar^ALP;

format longG

estar

rp

x

















