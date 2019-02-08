%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                           Two Equations Model
%                      Brianti, Cormun, February 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

% Structural Form System
% (1) Y_t = phi    *Y_{t-1} + psi     *X_t      + s_t^S
% (2) X_t = theta  *Y_{t-1} + (1-rho) *X_{t-1}  - s_t^Z

% Reduced Form Model
% (1) Y_t = (phi + psi*theta)  *Y_{t-1} + psi*(1-rho)   *X_{t-1} - psi*s_t^Z + s_t^S
% (2) X_t = theta              *Y_{t-1} + (1-rho)       *X_{t-1} -     s_t^Z

% Parameterization
phi   = 0.99;
psi   = -0.10;
rho   = 0.05;
thet  = .11;

% Eigenvalues Matrix Reduced Form
redu_matrix         = [phi+psi*thet psi*(1-rho); thet (1-rho)];
[AAA,eigens_matrix] = eig(redu_matrix);
eigens_matrix       = diag(eigens_matrix);
eigens_matrix       = sort(eigens_matrix)

% Parameterization
omeg = linspace(0,pi,10000);
omeg = omeg';
bet1 = phi + psi*thet + 1 - rho;
bet2 = - (phi)*(1-rho);
poly  = [1 -bet1 -bet2];
eigens_poly = sort(roots(poly));

if sum((eigens_matrix - eigens_poly).^2) > 10^(-12)
      error('Either Reduced Form Matrix or Polynomial in It is wrong')
end

DEN = 1 + bet1.^2 + bet2.^2 - 2.*bet1.*(1-bet2).*cos(omeg) - 2.*bet2.*cos(2*omeg);
NUM = 1 + (1-rho) - 2.*(1 - rho).*cos(omeg);
obj = NUM./DEN;
check = 2.*(1-rho).*sin(omeg).*DEN ...
      - NUM.*(4.*bet2.*sin(2.*omeg) + 2.*bet1.*(1 - bet2).*sin(omeg));

figure(1)
plot(obj)
legend('Objective')

figure(2)
plot(check)
legend('FOC')

% Sentiment Shocks
Ys(1) = 1;
Xs(1) = 0;
H = 35;
for i = 2:H
      Ys(i) = (phi + psi*thet)*Ys(i-1) + psi*(1-rho)*Xs(i-1);
      Xs(i) = thet*Ys(i-1) + (1-rho)*Xs(i-1);
end
% Technology Shocks
Yz(1) = -psi;
Xz(1) = -1;
for i = 2:H
      Yz(i) = (phi + psi*thet)*Yz(i-1) + psi*(1-rho)*Xz(i-1);
      Xz(i) = thet*Yz(i-1) + (1-rho)*Xz(i-1);
end
close all
figure('Position',[1 41 1920 963])
set(gcf,'color','w');
subplot(1,2,1)
plot(linspace(0,H-1,H),Ys,'--','linewidth', 2.5,'Color','r')
hold on
plot(linspace(0,H-1,H),Xs,'-','linewidth', 2,'Color','b')
grid on
plot([0:H-1]',0*[1:H]',':k','linewidth', 2);
xlabel('Period','fontsize',20);
ylabel('Level deviation from s.s.','fontsize',20);
title('Sentiment Shock', 'fontsize', 26);
hold off
subplot(1,2,2)
plot(linspace(0,H-1,H),Yz,'--','linewidth', 2.5,'Color','r')
hold on
plot(linspace(0,H-1,H),Xz,'-','linewidth', 2,'Color','b')
grid on
plot([0:H-1]',0*[1:H]',':k','linewidth', 2);
xlabel('Period','fontsize',20);
ylabel('Level deviation from s.s.','fontsize',20);
lgd = legend('Y','X','Location','South','Orientation','horizontal');
lgd.FontSize = 30;
legend('boxoff')
title('Technology Shock', 'fontsize', 26);
hold off

% Create the correct path
% base_path = pwd;
% warning off
% if exist([base_path '\Figures'], 'dir')
%       cd([base_path '\Figures']) %for Microsoft
% else
%       cd([base_path '/Figures']) %for Mac
% end
% if exist([base_path '\Export_Fig'], 'dir')
%       addpath([base_path '\Export_Fig']) %for Microsoft
% else
%       addpath([base_path '/Export_Fig']) %for Mac
% end
% warning on
% export_fig(['Theory_IRFs_BOTH.pdf'])
% cd(base_path) %back to the original path
% close


















