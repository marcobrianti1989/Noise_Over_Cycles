clear
close all
clc

% Find SS solution
fun   = @root2d;
x0    = [0.5,0.01]; % init values
x     = fsolve(fun,x0)

function F = root2d(X)

% Parameterization
del      = 0.05; % Capital depreciation
bet      = 0.98; % Discount factor
thetk    = 1;    % Technology of durables production
thetx    = 1;    % Technology of non-durables production
a        = 0.5;  % Weigh CES production of consumption good
v        = -1.5;  % Elasticity of substitution of CES production of consumption good
alp      = 0.6;  % Curvature production of non-durables
gam      = 0.97;  % Curvature production of durables
v0       = 1;    % Disutility of labor (linear)

% Key equations in ss

F(1) = (1 - bet*(1-del)) * v0 * ( thetk*gam*X(2)^(gam-1) )^(-1) ...
      - bet * ( a*( thetx*X(1)^alp )^v + (1-a)*( 1/del*thetk*X(2)^gam )^v )^(-1) ...
      * (1-a) * ( 1/del*thetk*X(2)^gam )^(v-1);


F(2) = v0 - ( a*( thetx*X(1)^alp )^v + (1-a)*( 1/del*thetk*X(2)^gam )^v )^(-1) ...
      * a * ( thetx*X(1)^alp )^(v-1) * thetx * alp * X(1)^(alp-1);

end






