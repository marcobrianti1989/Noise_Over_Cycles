clear 
close all

% From Beaudry, Galizia, Portier
% et = alp1*Xt + alp2*et-1 + alp3*et + mut

% Parameterization
alp1 = - 0.1;%19025; 
alp2 = .3;
del  = 0.01;
Xss  = 0;
ess  = 0;

bet1 = 1 - del + alp2;
bet2 = alp1 - alp2*(1 - del);         

matrix = [alp2 alp1; 1 (1-del)];
eig(matrix)

% Initialiation
H      = 50;
hor    = linspace(0,H,H);
e      = zeros(H,1);
SHOCK  = 1;

%-------------------------------------------------------------------------%
% IRFs to Sentiment Shock (Agents invest more today for no reasons)
%-------------------------------------------------------------------------%


e(1)   = SHOCK;
e(2)   = bet1*e(1);
for i = 3:H
      e(i)   = bet1*e(i-1) + bet2*e(i-2);
end

poly  = [1 -bet1 -bet2];
roots(poly)

figure('Position',[1 41 1920 963])
set(gcf,'color','w');
hold on
plot(hor,e,'linewidth',2,'color','r')
plot(hor,ones(1,H)*ess,'color','b')
title('Choice Variable','fontsize',20)
xlabel('Horizon','fontsize',16);
grid on
axis tight



asd


















