clear
close all

GAM2       = -1;
ETA2       = -1;
EPSI       = 0.1;
LAMBD      = 1;
ALP        = 2/3;    % Convexity of Production Function
DELT       = 0.05;   % Capital Depreciation Rate 
BET        = 0.99;   % Discount Factor
SIGM       = 0.5;    % IES

% Solution for Kss

check = ((1 - BET*R)./(GAM*LAMBD*R.^2)).^(1/ALP) - 5;