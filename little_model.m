clc
clear
close all

% param
del1   = 0.05;
del2   = 0.05;
mu1    = 0.1;
mu2    = -0.1;
mu3    = -0.5;
lam1   = 0.1;
lam2   = -0.2;
alp1   = 1.9;
alp2   = -0.2;

% System
A = 1 - del1 + lam1*alp1;
B = lam1*alp2 + lam2;
C = mu1*lam1*alp1 + mu2 + mu3*alp1;
D = 1 - del2 + mu1*lam1*alp2 + mu1*lam2 + mu3*alp2;
a = lam1;
b = 1;
c = mu1*lam1 + mu3;
d = mu1;

% Eigenvalues
MAT = [A B; C D]
eig(MAT)
