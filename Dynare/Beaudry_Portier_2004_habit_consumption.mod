% Brianti, Cormun; March 2019 - Beaudry, Portier (2004, JME)

%%%%%%% Defining Variables %%%%%%%

var 
C K LK LX X L I LOGTHETX LOGTHETK LOGSENT;          

%%%%% Aggregate Productivity Shock %%%%%

varexo eps_X eps_K eps_S;  

parameters
del bet thetk thetx a v alp gam gamL sig sigL v0 rhox rhok rhos extX extK; 

% Parameterization
del      = 0.05;  % Capital depreciation
bet      = 0.98;  % Discount factor
thetk    = 1;     % Technology of durables production
thetx    = 1;     % Technology of non-durables production
a        = 0.5;   % Weigh CES production of consumption good
v        = -1.5;  % Elasticity of substitution of CES production of consumption good
alp      = 0.6;   % Curvature production of non-durables
gam      = 0.97;  % Curvature production of durables
gamL     = 1;     % frish elasticity
sig      = 0;     % degree of consumption habit
sigL     = 0.999; % degree of leisure habit
v0       = 1;     % Disutility of labor (linear)
rhox     = 0.5;   % Persistence of tech shock of durables production
rhok     = 0;     % Persistence of tech shock of non-durables production
rhos     = 0.5;   % Persistence of sentiment shock
extX     = 0;     % externality in production of X
extK     = 0;     % externality in production of K

% Defining functionals forms and derivatives

model; 

exp(LOGSENT) * v0/(1 - L + sigL*L(-1))^gamL * ( exp(LOGTHETK) * LK^extK * thetk * gam * LK^(gam-1) )^(-1) = bet * ( (1-del) * exp(LOGSENT(+1)) * v0/(1 - L(+1) + sigL*L)^gamL * ( exp(LOGTHETK(+1))*LK(+1)^extK*thetk*gam*LK(+1)^(gam-1) )^(-1) + C(+1)^(1-v)/(C(+1) - sig*C) * (1-a) * K(+1)^(v-1) );

exp(LOGSENT) * v0/(1 - L + sigL*L(-1))^gamL =  1/(C - sig*C(-1)) * C^(1-v) * a * ( exp(LOGTHETX) * LX^extX * thetx * LX^alp )^(v-1) * exp(LOGTHETX) * LX^extX * thetx * alp * LX^(alp-1);

C = ( a*( X )^v + (1-a)*( K )^v )^(1/v);

L = LK + LX;

K = (1-del)*K(-1) + I(-1);

I = exp(LOGTHETK) * LK^extK * thetk * LK^gam; 

X = exp(LOGTHETX)* LX^extX * thetx * LX^alp;

LOGTHETX  = rhox*LOGTHETX(-1) + eps_X;

LOGTHETK = rhok*LOGTHETK(-1) + eps_K;

LOGSENT = rhos*LOGSENT(-1) - eps_S;

end;

LOGTHETXss = 0;
LOGTHETKss = 0;
LOGSENTss  = 0;
LXss       = 0.51199;
LKss       = 0.10104;
Xss        = thetx * LXss^alp;
Kss        = 1/del * thetk * LKss^gam;
Css        = ( a*( Xss )^v + (1-a)*( Kss )^v )^(1/v);
Lss        = LKss + LXss;
Iss        = thetk * LKss^gam; 

                           
%%%%% Initialization %%%%%
initval;

LOGTHETX   = LOGTHETXss;
LOGTHETK   = LOGTHETKss;
LOGSENT    = LOGSENTss;
LX         = LXss;
LK         = LKss;
X          = Xss;
K          = Kss;
C          = Css;
L          = Lss;
I          = Iss; 

end;


shocks;
  var eps_X    = 1;
  %var eps_K    = 1;
  var eps_S    = 1;
end;

steady;
check;

stoch_simul(irf=24, order=1) LOGSENT LOGTHETX C I L X LK LX K;

%stoch_simul(periods=100000, hp_filter = 1600, order=2,ar = 0, nofunctions,nograph,nodecomposition,nocorr) jobphim logz logR logU logV logFPHIC;



