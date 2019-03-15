% Brianti, Cormun January 2019 - Buaudry, Galizia, Portier (2019, AER)

%%%%%%% Defining Variables %%%%%%%

var 
                              
Z                      % (1) Endogenous Technology
K                      % (2) Level of Capital
R                      % (3) Default Rate
EPS                    % (4) Weights for new Tech level          
I                      % (5) Investment
Y                      % (6) Output
Rf                     % (7) Risk free
RHO                    % (8) Probability of default

%%%%% Aggregate Productivity Shock %%%%%

varexo eps;  

% Parameterization
parameters
BET GAM ALP DELT ETA1 ETA2 EPSI LAMBD SIGM;

GAM1       = 1;
GAM2       = -1;
ETA2       = -1;
EPSI       = 0.1;
LAMBD      = 1;
ALP        = 2/3;    % Convexity of Production Function
DELT       = 0.05;   % Capital Depreciation Rate 
BET        = 0.99;   % Discount Factor
SIGM       = 0.5;    % IES

% Defining functional forms and derivatives

model; 

1 = BET*( ( Y - I + Rf(-1)*I(-1) ) / ( Y(+1) - I(+1) + Rf*I ) )^SIGM *Rf;   %(1)

Y = Z * K^ALP;  %(2)

I = ETA1*K(-1) + ETA2*R; %(3)

( 1 - RHO ) * R = Rf;  %(4)

K = ( 1 - DEL ) * K(-1) + ( 1 - RHO ) * I;         %(5)

RHO = GAM * Y;   %(6)

Z = ( 1 - EPS ) * Z(-1) + EPS * LAMBD * R + epsZ;    %(7)

EPS = ( 1 - RHO ) * I / K;

R = R + epsS;

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


stoch_simul(irf=30, order=1) logthet logzeta e x c;

stoch_simul(irf=20, order=1) logthet logzeta e x rp;


%stoch_simul(periods=100000, hp_filter = 1600, order=2,ar = 0, nofunctions,nograph,nodecomposition,nocorr) jobphim logz logR logU logV logFPHIC;



