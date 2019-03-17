clear
close all



% Solution for Kss



fun = @root2d;
x0 = [0.1,0.1,0.1];
x = fsolve(fun,x0);


function F = root2d(X)

Z          = 2;      % mean of productivity of all the firms 
XX         = 1;   % fixed cost of production
EPS        = 0.9;    % percentage of capital recover if default
FF         = 0.1;      % fixed cost of recovering capital if default
ALP        = 2/3;    % Convexity of Production Function
BET        = 0.99;   % Discount Factor
GAM        = 0.1;    % habit consumption formation
DEL        = 0.05;   % Capital depreciation 
Rss        = 1/BET - 1 + DEL;

F(1) = Rss.*X(2) - EPS*Z*X(1).*X(2).^ALP + FF;
F(2) = Z*X(1).*X(2).^ALP - Rss.*X(2) - XX;
F(3) = 1/2*Z*(1 - X(1).^2).*X(2).^ALP + (1 - X(1)).*( (1-DEL)*X(2) - XX ) ...
      - X(3) - (1-X(1)).*X(2);

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