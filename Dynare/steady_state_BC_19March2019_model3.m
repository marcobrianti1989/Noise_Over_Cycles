clear
close all



% Solution ss
fun = @root2d;
x0 = [0.1,0.1];
x = fsolve(fun,x0)


function F = root2d(X)

Z          = 1;      % mean of productivity of all the firms 
ETA        = 0.9;    % elasticity of financial fragility to past output
EPS        = 0.01;    % percentage of capital recover if default
FF         = 1;    % fixed cost of recovering capital if default
ALP        = 2/3;    % Convexity of Production Function
BET        = 0.99;   % Discount Factor
GAM        = 0.5;    % habit consumption formation
DEL        = 0.05;   % Capital depreciation 
SIG        = 2;      % IES
PSI        = 1;      % multiplicative disutility of labor
CHI        = 2;      % Frish elasticity

Rss        = 1/BET - 1 + DEL; % R in ss

F(1) = PSI*X(2)^(1+CHI)*(1/(1-ALP)) ...
      - EPS*Z*(ALP/(1-ALP)*PSI*X(2)^(1+CHI))^ALP*X(2)^(1-ALP) + FF;
F(2) = X(1) + DEL/Rss*ALP/(1-ALP)*PSI*X(2)^(1+CHI) ...
      - Z*(ALP/(1-ALP)*PSI*X(2)^(1+CHI))^ALP*X(2)^(1-ALP);

end

% rho = linspace(0,1,100);
% n   = 0.75;
% 
% obj1 = PSI*(1 - rho.^2).*n.^(1+CHI) - EPS*Z*rho.*n.^ALP + F;
% 
% obj2 = rho.*n.^ALP - PSI*(1 - rho.^2).*n.^(1+CHI) - X./(1/2*Z*(1-rho.^2).*n.^ALP);
% 
% plot(obj1)
% hold on
% plot(obj2)