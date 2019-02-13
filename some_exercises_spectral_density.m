% omeg = linspace(0,2*pi,5000);
%
% Z = exp(1i.*omeg)'
% Z2 = cos(omeg) + 1i*sin(omeg);
% plot(Z)
% hold on
% plot(Z2)

clear
close all

T = 1000000;
alp1 = 0.9;
alp2 = -0.80;
alp3 = 0;
alp4 = 0;
sig  = 0.01;
X(1,1) = sig*rand(1,1);
X(2,1) = alp1*X(1,1) + sig*randn(1,1);
X(3,1) = alp1*X(2,1) + alp2*X(1,1) + sig*randn(1,1);
X(4,1) = alp1*X(3,1) + alp2*X(2,1) + alp3*X(1,1) + sig*randn(1,1);
for ii = 5:T
      X(ii,1) = alp1*X(ii-1,1) + alp2*X(ii-2,1) + alp3*X(ii-3,1) + alp4*X(ii-4,1) + sig*randn(1,1);
end
% plot(X)
% close
omeg = linspace(0,pi,50)';
Z    = exp(-1i*omeg);
lags = 100;
for iZ = 1:length(omeg)
      [GxZ(iZ,1), Cx] = autocov_GenFunc(X,lags,Z(iZ));
end
periodicity = (2*pi./omeg);
figure(1)
plot(omeg,GxZ)
figure(2)
plot(periodicity,GxZ)

