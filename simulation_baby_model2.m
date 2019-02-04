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

matrix = [alp2 alp1; 1 (1-del)];
eig(matrix)

% Initialiation
H      = 50;
hor    = linspace(0,H,H);
e      = zeros(H,1);
X      = zeros(H,1);
SHOCK  = 1;
X(1)   = Xss; % state variable always Xss on impact

%-------------------------------------------------------------------------%
% IRFs to Sentiment Shock (Agents invest more today for no reasons)
%-------------------------------------------------------------------------%


e(1)   = SHOCK;
X(1)   = Xss;
for i = 2:H
      % Choice Variable
      e(i)   = alp2*e(i-1) + alp1*X(i-1);
      % LOM of X
      X(i)   = e(i-1) + (1 - del)*X(i-1);    
end

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

clear 
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  Spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameterization
omeg = linspace(0,pi,10000);
omeg = omeg';
phi = 0.9;
psi = -0.11;
del = 0.05;
bet1 = phi + psi + 1 - del;
bet2 = phi*(1-del);

DEN = 1 + bet1.^2 + bet2.^2 - 2.*bet1.*(1-bet2).*cos(omeg) - 2.*bet2.*cos(2*omeg);
NUM = 1 + (1-del) - 2.*(1 - del).*cos(omeg);
obj = NUM./DEN;
check = 2.*(1-del).*sin(omeg).*DEN ...
      - NUM.*(4.*bet2.*sin(2.*omeg) + 2.*bet1.*(1 - bet2).*sin(omeg));

figure(1)
plot(obj)
legend('Objective')

figure(2)
plot(check)
legend('FOC')

























