clear 
close all

% From Beaudry, Galizia, Portier
% et = alp0 + alp1*Xt + alp2*et-1 + alp3*et + mut

% Parameterization
alp0 = 0.4;
alp1 = - 0.1; 
alp2 = 1;
del  = 0.05;

% Steady State
eold  = 0;
Xold  = 0;
err   = 1;
simss = 0;
while err > 10^(-14) && simss < 1000
      % Choice Variable
      enew   = alp0 + alp1*Xold + alp2*eold;
      % LOM of X
      Xnew   = (1 - del)*Xold + enew;    
      err    = abs(Xnew - Xold) + abs(enew - eold);
      simss  = simss + 1;
      eold   = enew;
      Xold   = Xnew;
end
Xss  = Xnew;
ess  = enew;

% Initialiation
H      = 20;
hor    = linspace(0,H,H);
e      = zeros(H,1);
X      = zeros(H,1);
SHOCK  = 1;
X(1)   = Xss; % state variable always Xss on impact

%-------------------------------------------------------------------------%
% IRFs to Sentiment Shock (Agents invest more today for no reasons)
%-------------------------------------------------------------------------%


e(1)   = alp0 + alp1*Xss + alp2*ess + SHOCK;
X(2)   = (1 - del)*Xss + e(1);
for i = 2:H
      % Choice Variable
      e(i)   = alp0 + alp1*X(i) + alp2*e(i-1);
      % LOM of X
      X(i+1) = (1 - del)*X(i) + e(i);    
end
X = X(1:end-1);

figure('Position',[1 41 1920 963])
set(gcf,'color','w');
subplot(1,2,1)
hold on
plot(hor,e,'linewidth',2,'color','r')
plot(hor,ones(1,H)*ess,'color','b')
title('Choice Variable','fontsize',20)
xlabel('Horizon','fontsize',16);
grid on
axis tight
hold off
subplot(1,2,2)
hold on
plot(hor,X,'linewidth',2,'color','r')
plot(hor,ones(1,H)*Xss,'color','b')
title('State Variable','fontsize',20)
xlabel('Horizon','fontsize',16);
grid on
axis tight
hold off

%-------------------------------------------------------------------------%
% IRFs to Fundamental Shock (Exogenous increase in capital today)
%-------------------------------------------------------------------------%

e(1)   = alp0 + (alp1)*Xss + alp2*ess;
X(2)   = (1 - del)*Xss + e(1) + SHOCK;
for i = 2:H
      % Choice Variable
      e(i)   = alp0 + (alp1)*X(i) + alp2*e(i-1);
      % LOM of X
      X(i+1) = (1 - del)*X(i) + e(i);    
end
X = X(1:end-1);

figure('Position',[1 41 1920 963])
set(gcf,'color','w');
subplot(1,2,1)
hold on
plot(hor,e,'linewidth',2,'color','r')
plot(hor,ones(1,H)*ess,'color','b')
title('Choice Variable','fontsize',20)
xlabel('Horizon','fontsize',16);
grid on
axis tight
hold off
subplot(1,2,2)
hold on
plot(hor,X,'linewidth',2,'color','r')
plot(hor,ones(1,H)*Xss,'color','b')
title('State Variable','fontsize',20)
xlabel('Horizon','fontsize',16);
grid on
axis tight
hold off

