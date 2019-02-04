clear
close all

% Structural Form System
% (1) I_t = phi*I_{t-1} + psi*K_t + s_t^I
% (2) K_t = I_{t-1} + (1-del)*K_{t-1} + s_t^K

% Reduced Form Model
% (1) I_t = (phi + psi)*I_{t-1} + psi*(1-del)*K_t + psi*s_t^K + s_t^I
% (2) K_t = I_{t-1} + (1-del)*K_{t-1} + s_t^K

% Parameterization
phi = 0.4;
psi = -0.1;
del = 0.05;

% Eigenvalues Matrix Reduced Form
redu_matrix       = [phi+psi psi*(1-del); 1 (1-del)];
[AAA,eigens_matrix] = eig(redu_matrix);
eigens_matrix = diag(eigens_matrix);
eigens_matrix = sort(eigens_matrix);

% Parameterization
omeg = linspace(0,pi,10000);
omeg = omeg';
bet1 = phi + psi + 1 - del;
bet2 = - (phi)*(1-del);
poly  = [1 -bet1 -bet2];
eigens_poly = sort(roots(poly));

if sum((eigens_matrix - eigens_poly).^2) > 10^(-12)
      error('Either Reduced Form Matrix or Polynomial in It is wrong')
end

DEN = 1 + bet1.^2 + bet2.^2 - 2.*bet1.*(1-bet2).*cos(omeg) - 2.*bet2.*cos(2*omeg);
NUM = 1 + (1-del) - 2.*(1 - del).*cos(omeg);
obj = NUM./DEN;
check = 2.*(1-del).*sin(omeg).*DEN ...
      - NUM.*(4.*bet2.*sin(2.*omeg) + 2.*bet1.*(1 - bet2).*sin(omeg));

figure(1)
plot(obj)
legend('Objective')
asd
figure(2)
plot(check)
legend('FOC')