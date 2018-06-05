%-----------------
% Housekeeping    
%----------------

clear all;  
close all;
clc; 

addpath('./auxfiles')
addpath('./results')

%--------------------
% Model Specification  
%--------------------
do_Figure       = 1;
i_var_instr     = {'TAXNARRATIVE'};                           % Instruments for proxy SVAR
VAR.select_vars = {'TAX_S','G_S','GDP_S','TB3MS','CPI_PIQ4'}; % Variables in the VAR
p               = 4; % Number of lags
nex_            = 1; % Constant
n               = size(VAR.select_vars,2);                    % Number of variables
nP              = 2; % Number of fiscal variables
ivar_t          = 1; % Position of tax revenue
ivar_g          = 2; % Position of government spending
T0              = p; % Size of pre-sample for Minnesota Prior
Horizon         = 1; % Horizon for Impulse Responses (Can only be set to 1 in this replication code)

%--------------------
% Load Data  
%--------------------

load DATASET
GDPRATIO = DATASET.GDPRATIO;

YY           = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,VAR.select_vars)));
VAR.ivar_str = DATASET.FIGLABELS(cell2mat(values(DATASET.MAP,VAR.select_vars)),1);
VAR.MAP      = containers.Map([VAR.select_vars],1:size(YY,2));
VAR.mm       = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,i_var_instr)));
M            = VAR.mm(T0+1:end,:);

%-------------------------
% Generate Minnesota Prior
%-------------------------

vm_dummy;

%-------------------------
% Estimation Preliminaries
%-------------------------

J = [eye(n);repmat(zeros(n),p-1,1)]; 
F = zeros(n*p,n*p);  % Matrix for Companion Form
I  = eye(n);

for i=1:p-1
    F(i*n+1:(i+1)*n,(i-1)*n+1:i*n) = I;
end

X = [XXdum; XXact];
Y = [YYdum; YYact];

T = size(X, 1);
ndum = size(XXdum, 1);

% OLS estimates

B = (X'*X)\(X'*Y); % Point estimates
U = Y-X*B;      % Residuals
Sigmau = U'*U/(T-p*n-1);   % Covariance matrix of residuals
F(1:n,1:n*p)    = B(1:n*p,:)';
sTT = Sigmau(1,1); sYY = Sigmau(3,3); sGG = Sigmau(2,2);
sYT = Sigmau(1,3); sYG = Sigmau(2,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define objects to store results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select grid for elasticities of taxes (etaTY) and spending (etaGY)

etaTYgrid = -2:0.01:5; 
etaGYgrid = -2.5:0.01:4.5;

nGrid = length(etaTYgrid);

F21TOrd1 = zeros(nGrid,1);
F21GOrd1 = zeros(nGrid,1);
eeta = zeros(n-2,nP);

for ii = 1:nGrid

    eeta(1,1)=etaTYgrid(ii);
    eeta(1,2)=etaGYgrid(ii);

    % Compute tax and spending multipliers when shocks are ordered first
    F21Ord1 = vm_mult_simple(F,J,Horizon,eeta,U,T,n,p,nP,DATASET.GDPRATIO);
    F21TOrd1(ii,1) = F21Ord1(1);
    F21GOrd1(ii,1) = F21Ord1(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute elasticities and multipliers to plot in Figure 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------
%Blanchard and Perotti
%---------------------

etaTYBP = 1.7; etaGYBP = 0;
F21TBP = -(sYT - etaTYBP*sYY)/(etaTYBP^2*sYY + sTT -2*etaTYBP*sYT)/GDPRATIO(2);
F21GBP =  (sYG - etaGYBP*sYY)/(etaGYBP^2*sYY + sGG -2*etaGYBP*sYG)/GDPRATIO(3);

BPT = [etaTYBP F21TBP];
BPG = [etaGYBP F21GBP];

%-------------------------
%Penalty Function Approach
%-------------------------
nperPF = 1;
[ffactorPF, A0PF, IRFPF] = penalty_fun_simple(Sigmau,n,p,F,J,Horizon,nperPF);
% Store output elasticities of taxes and spending;
etaPF_T = -A0PF(3,3)/A0PF(1,3);
etaPF_G = -A0PF(3,2)/A0PF(2,2);
% Impact Tax and Spending Multipliers from Analytical Solution;
F21TPF = -(sYT - etaPF_T*sYY)/(etaPF_T^2*sYY + sTT -2*etaPF_T*sYT)/GDPRATIO(2);
F21GPF =  (sYG - etaPF_G*sYY)/(etaPF_G^2*sYY + sGG -2*etaPF_G*sYG)/GDPRATIO(3);

PFAT = [etaPF_T F21TPF];
PFAG = [etaPF_G F21GPF];

%-------------------
%Proxy SVAR Approach
%-------------------

F1 = (M'*M)\(M'*U(ndum+1:end,:));
xi = F1/F1(1);
% Calculate output elasticity of taxes;
etaPROXY = (sYT-xi(3)*sTT)/(sYY-xi(3)*sYT);
F21TPROXY = -(sYT - etaPROXY*sYY)/(etaPROXY^2*sYY + sTT -2*etaPROXY*sYT)/GDPRATIO(2);

PROXY = [etaPROXY F21TPROXY];

%-------------------
%Cholesky
%-------------------

etaCholesky_T = sYT/sYY;
etaCholesky_G = sYG/sYY;
CHOLESKYG = [etaCholesky_G 0];
CHOLESKYT = [etaCholesky_T 0];

%-------------------
% Produce Figure
%-------------------

fig = figure(1);
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
a = 56;

subplot(121)
axis([-2 5 -2.2 2.2])
vline(0,'k');
hold on
hline(0,'k');
plot(etaTYgrid,F21TOrd1(:,1),'k','linewidth',1.5);
h1 = scatter(CHOLESKYT(1),CHOLESKYT(2),a*2,'p','filled','MarkerEdgeColor','r','MarkerFaceColor','r');
h2 = scatter(BPT(1),BPT(2),a,'s','filled','MarkerEdgeColor','b','MarkerFaceColor','b');
h3 = scatter(PFAT(1),PFAT(2),a,'o','MarkerEdgeColor','k','LineWidth',2);
h4 = scatter(PROXY(1),PROXY(2),a,'d','filled','MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor',[0 0.6 0]);
xlab = sprintf('Tax Rule Coefficient (\\psi^{tr}_{gdp})');
xlabel(xlab,'FontSize',12)
ylabel('Impact Tax Multiplier','FontSize',12)
set(gca,'YTick',-2:1:2,'YTickLabel',{'-2.0','-1.0','0.0','1.0','2.0'})
set(gca,'box','on')

subplot(122)
axis([-2.5 2.5 -2.2 2.2])
vline(0,'k');
hold on
hline(0,'k');
plot(etaGYgrid,F21GOrd1(:,1),'k','linewidth',1.5);
scatter(CHOLESKYG(1),CHOLESKYG(2),a*2,'p','filled','MarkerEdgeColor','r','MarkerFaceColor','r');
scatter(BPG(1),BPG(2),a,'s','filled','MarkerEdgeColor','b','MarkerFaceColor','b');
scatter(PFAG(1),PFAG(2),a,'o','MarkerEdgeColor','k','LineWidth',1.5);
xlab = sprintf('Spending Rule Coefficient (\\psi^{g}_{gdp})');
xlabel(xlab,'FontSize',12)
ylabel('Impact Spending Multiplier','FontSize',12)
set(gca,'YTick',-2:1:2,'YTickLabel',{'-2.0','-1.0','0.0','1.0','2.0'})
set(gca,'box','on')


lab1 = sprintf('Cholesky ($\\psi^{p,chol}_{gdp}$)');
lab2 = sprintf('Blanchard-Perotti{   }');
lab3 = sprintf('Penalty Function{   }');
lab4 = sprintf('Proxy SVAR{   }');
legend([h1,h2,h3,h4],{lab1,lab2,lab3,lab4},'Orientation','Horizontal','Location',[0.35,0,0.3,0.05],'FontSize',10,'box','off','Interpreter','latex')  

% Print figure

dim = [18,9.5];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(fig,'-depsc','./results/Figure1');
print(fig,'-dpdf' ,'./results/Figure1');
