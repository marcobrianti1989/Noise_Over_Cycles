clear
close all



% Solution for Kss



fun = @root2d;
x0 = [0.1,0.1];
x = fsolve(fun,x0);


function F = root2d(X)

Z          = 2;      % mean of productivity of all the firms 
XX         = 0.01;    % fixed cost of production
EPS        = 1;      % percentage of capital recover if default
FF         = 0.5;    % fixed cost of recovering capital if default
ALP        = 2/3;    % Convexity of Production Function
BET        = 0.99;   % Discount Factor
PSI        = 0.1;    % linear consumption disutility
CHI        = 2;      % Frish

F(1) = PSI*(1 - X(1).^2).*X(2).^(1+CHI) - EPS*Z*X(1).*X(2).^ALP + FF;
F(2) = X(1).*X(2).^ALP - PSI*(1 - X(1).^2).*X(2).^(1+CHI) - XX./(1/2*Z*(1-X(1).^2).*X(2).^ALP);

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