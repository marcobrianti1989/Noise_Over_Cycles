tic

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
do_Figure       = 2;
i_var_instr     = {'TAXNARRATIVE'};           % Instrument for proxy SVAR
VAR.select_vars = {'TAX_S','G_S','GDP_S','CPI_PIQ4','TB3MS'};
p               = 4;      % Number of lags
nex_            = 1;      % Constant
n               = size(VAR.select_vars,2); % Number of variables
nP              = 2; % Number of fiscal variables
T0              = p;      % Size of pre-sample for Minnesota Prior
Horizon         = 41;     % Horizon for Impulse Responses
nd              = 20000;  % Number of draws in MC chain
bburn           = 0.2*nd; % Share of draws to burn

ptileVEC = [0.16 0.50 0.84]; % Percentiles for posterior distributions

do_BP        = 1;   % Execute Blanchard Perotti (2002) identification
do_PF        = 1;   % Execute penalty function identification
do_proxySVAR = 1;   % Execute proxy SVAR identification

% Fiscal elasticities for Blanchard Perotti identification
BP_EETA_SIMPLE      = zeros(n-2,nP);
BP_EETA_SIMPLE(1,1) = 1.7;
BP_EETA_GENERAL     = [1.7 0; 1.25 -0.50; 0 0];

% Number of periods for penalty function identification
nperPF = 1;


%--------------------
% Load Data  
%--------------------

load DATASET
GDPRATIO = DATASET.GDPRATIO;

YY                   = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,VAR.select_vars)));
VAR.ivar_str         = DATASET.FIGLABELS(cell2mat(values(DATASET.MAP,VAR.select_vars)),1);
VAR.MAP              = containers.Map([VAR.select_vars],[1:size(YY,2)]);
VAR.mm               = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,i_var_instr)));
M                    = VAR.mm(T0+1:end,:);
nIV                  = size(VAR.mm,2);


%-------------------------
% Generate Minnesota Prior
%-------------------------

vm_dummy;

%-------------------------
% Estimation Preliminaries
%-------------------------

J = [eye(n);repmat(zeros(n),p-1,1)]; % Page 12 RWZ
F = zeros(n*p,n*p);    % Matrix for Companion Form
I  = eye(n);
for i=1:p-1
    F(i*n+1:(i+1)*n,(i-1)*n+1:i*n) = I;
end

X = [XXdum; XXact];
Y = [YYdum; YYact];
T = size(X, 1);
ndum = size(XXdum, 1);

% Compute OLS estimates
B = (X'*X)\(X'*Y);       % Point estimates
U = Y-X*B;               % Residuals
Sigmau = U'*U/(T-p*n-1); % Covariance matrix of residuals
F(1:n,1:n*p)    = B(1:n*p,:)';

%------------------------------------------------------------
% Set preliminaries for priors
%------------------------------------------------------------

N0=zeros(size(X',1),size(X,2));
nnu0=0;
nnuT = T +nnu0;
NT = N0 + X'*X;    
Bbar0=B;
S0=Sigmau;
BbarT = NT\(N0*Bbar0 + (X'*X)*B);
ST = (nnu0/nnuT)*S0 + (T/nnuT)*Sigmau + (1/nnuT)*((B-Bbar0)')*N0*(NT\eye(n*p+nex_))*(X'*X)*(B-Bbar0); %% Constant (check)
STinv = ST\eye(n);

%------------------------------------------------------------
% Define objects that store the draws
%------------------------------------------------------------
BP_DM_T_SIMPLE = zeros(nd-bburn,Horizon);
BP_DM_G_SIMPLE = zeros(nd-bburn,Horizon);
BP_DM_T_GENERAL = zeros(nd-bburn,Horizon);
BP_DM_G_GENERAL = zeros(nd-bburn,Horizon);
BP_EETA_TG = zeros(nd-bburn,1);
PF_DM_T_SIMPLE = zeros(nd-bburn,Horizon);
PF_DM_G_SIMPLE = zeros(nd-bburn,Horizon);
PF_DM_T_GENERAL = zeros(nd-bburn,Horizon);
PF_DM_G_GENERAL = zeros(nd-bburn,Horizon);
PF_EETAT_GENERAL = zeros(nd-bburn,n,1);
PF_EETAG_GENERAL = zeros(nd-bburn,n,1);
PROXY_DM_T_SIMPLE = zeros(nd-bburn,Horizon);
PROXY_DM_T_GENERAL = zeros(nd-bburn,Horizon);
PROXY_EETA_GENERAL = zeros(nd-bburn,n-1,1);

PF_EETA_SIMPLE      = zeros(n-2,nP);
PROXY_EETA_SIMPLE   = zeros(n-2,nP);
record=0;     
counter = 0;


disp('                                                                  ');
disp('        BAYESIAN ESTIMATION OF VAR VIA BLOCK MCMC                 ');
disp('                                                                  ');

%-----------------
% MCMC Chain 
%----------------


while record<nd

    %------------------------------------------------------------
    % Gibbs Sampling Algorithm
    %------------------------------------------------------------
    % STEP ONE: Draw from the B, SigmaB | Y
    %------------------------------------------------------------
    % Step 1: Draw from the marginal posterior for Sigmau p(Sigmau|Y,X)

    R=mvnrnd(zeros(n,1),STinv/nnuT,nnuT)';
    Sigmadraw=(R*R')\eye(n);

    % Step 2: Taking newSigma as given draw for B using a multivariate normal    
    bbeta = B(:);
    SigmaB = kron(Sigmadraw,NT\eye(n*p+nex_));
    SigmaB = (SigmaB+SigmaB')/2;
    Bdraw = mvnrnd(bbeta,SigmaB);
    % Storing unrestricted draws
        
    Bdraw = reshape(Bdraw,n*p+nex_,n);% Reshape Bdraw from vector to matrix
    Udraw = Y-X*Bdraw;      % Store residuals for IV regressions
    
    F(1:n,1:n*p)    = Bdraw(1:n*p,:)';

    record=record+1;
    counter = counter +1;
    if counter==0.05*nd
        disp(['         DRAW NUMBER:   ', num2str(record)]);
        disp('                                                                  ');
        disp(['     REMAINING DRAWS:   ', num2str(nd-record)]);
        disp('                                                                  ');
        counter = 0;
    end

    if record > bburn
        
        if do_BP == 1
            F21Ord1BP = vm_mult_simple(F,J,Horizon,BP_EETA_SIMPLE,Udraw,T,n,p,nP,GDPRATIO);
            BP_DM_T_SIMPLE(record-bburn,:) = F21Ord1BP(:,1);
            BP_DM_G_SIMPLE(record-bburn,:) = F21Ord1BP(:,2);

            [F21Ord1BPfull,BP_EETA_TG(record-bburn,:)] = vm_mult_bp(F,J,Horizon,BP_EETA_GENERAL,Udraw,T,n,p,GDPRATIO);
            BP_DM_G_GENERAL(record-bburn,:) = F21Ord1BPfull(:,1);
            BP_DM_T_GENERAL(record-bburn,:) = F21Ord1BPfull(:,2);
        end
        
        if do_PF == 1
            [~, A0PF, ~] = penalty_fun_simple(Sigmadraw,n,p,F,J,Horizon,nperPF);
            PF_EETA_SIMPLE(1,1) = -A0PF(3,3)/A0PF(1,3); % First row tax, third row GDP; Second column is tax shock;
            PF_EETA_SIMPLE(1,2) = -A0PF(3,2)/A0PF(2,2); % Second row spending, third row GDP; Third column is spending shock;
            F21Ord1PF = vm_mult_simple(F,J,Horizon,PF_EETA_SIMPLE,Udraw,T,n,p,nP,GDPRATIO);
            PF_DM_T_SIMPLE(record-bburn,:) = F21Ord1PF(:,1);
            PF_DM_G_SIMPLE(record-bburn,:) = F21Ord1PF(:,2);
            
            
            [~, A0PF, IRFPF] = penalty_fun_general(Sigmadraw,n,p,F,J,Horizon,nperPF);
            PF_EETAT_GENERAL(record-bburn,:) = -A0PF(:,5)/A0PF(1,5);
            PF_EETAG_GENERAL(record-bburn,:) = -A0PF(:,4)/A0PF(2,4);
            PF_DM_T_GENERAL(record-bburn,:) = -squeeze(IRFPF(:,3,5))*A0PF(1,5)/GDPRATIO(2);
            PF_DM_G_GENERAL(record-bburn,:) =  squeeze(IRFPF(:,3,4))*A0PF(2,4)/GDPRATIO(3);           
        end
                
        if do_proxySVAR == 1
            sTT = Sigmadraw(1,1); sYY = Sigmadraw(3,3); sYT = Sigmadraw(1,3);
            F1 = (M'*M)\(M'*Udraw(ndum+1:end,:));
            xi = F1/F1(1);
            PROXY_EETA_SIMPLE(1,1) = (sYT-xi(3)*sTT)/(sYY-xi(3)*sYT);
            F21Ord1 = vm_mult_simple(F,J,Horizon,PROXY_EETA_SIMPLE,Udraw,T,n,p,nP,GDPRATIO);
            PROXY_DM_T_SIMPLE(record-bburn,:) = F21Ord1(:,1);

            [PROXY_EETA, PROXY_IRF_T_GENERAL] =  proxy_svar_fiscal(Sigmadraw,M,Udraw(ndum+1:end,:),n,p,F,J,Horizon);
            PROXY_EETA_GENERAL(record-bburn,:) = PROXY_EETA';
            PROXY_DM_T_GENERAL(record-bburn,:) = PROXY_IRF_T_GENERAL/GDPRATIO(2);

        end
    end
end

SVAR.BP_DM_T_SIMPLE_Q     = quantile(BP_DM_T_SIMPLE,ptileVEC);
SVAR.BP_DM_G_SIMPLE_Q     = quantile(BP_DM_G_SIMPLE,ptileVEC);
SVAR.BP_DM_T_GENERAL_Q    = quantile(BP_DM_T_GENERAL,ptileVEC);
SVAR.BP_DM_G_GENERAL_Q    = quantile(BP_DM_G_GENERAL,ptileVEC);
SVAR.BP_EETA_TG_Q         = quantile(BP_EETA_TG,ptileVEC);
SVAR.PF_DM_T_SIMPLE_Q     = quantile(PF_DM_T_SIMPLE,ptileVEC);
SVAR.PF_DM_G_SIMPLE_Q     = quantile(PF_DM_G_SIMPLE,ptileVEC);
SVAR.PF_DM_T_GENERAL_Q    = quantile(PF_DM_T_GENERAL,ptileVEC);
SVAR.PF_DM_G_GENERAL_Q    = quantile(PF_DM_G_GENERAL,ptileVEC);
SVAR.PF_EETAT_GENERAL_Q   = quantile(PF_EETAT_GENERAL,ptileVEC);
SVAR.PF_EETAG_GENERAL_Q   = quantile(PF_EETAG_GENERAL,ptileVEC);
SVAR.PROXY_DM_T_SIMPLE_Q  = quantile(PROXY_DM_T_SIMPLE,ptileVEC);
SVAR.PROXY_DM_T_GENERAL_Q = quantile(PROXY_DM_T_GENERAL,ptileVEC);
SVAR.PROXY_EETA_GENERAL_Q = quantile(PROXY_EETA_GENERAL,ptileVEC);


savefileirf = strcat('./results/Results_Figure2Table1','.mat');
save(savefileirf,'SVAR'); 

toc

%% Figure 2

transp90 = 0.25;
font_num = 12;
linW = 2;

fig = figure(1);
subplot(321)
hline(0,'k');
hold on
a = SVAR.BP_DM_T_GENERAL_Q(1,:);
b = SVAR.BP_DM_T_GENERAL_Q(3,:);
[~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
hold on
h1 = plot(0:1:Horizon-1,SVAR.BP_DM_T_GENERAL_Q(2,:),'b','LineWidth',linW);
h2 = plot(0:1:Horizon-1,SVAR.BP_DM_T_SIMPLE_Q(2,:),'r--','LineWidth',linW);
title('Blanchard-Perotti: Tax Multiplier','FontSize',font_num)
set(gca,'YTick',-0.5:0.5:2.5)
set(gca,'XTick',0:4:Horizon-1)
axis([0 40 -0.75 2.5])
box on

subplot(322)
hline(0,'k');
hold on
a = SVAR.BP_DM_G_GENERAL_Q(1,:);
b = SVAR.BP_DM_G_GENERAL_Q(3,:);
[~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
hold on
plot(0:1:Horizon-1,SVAR.BP_DM_G_GENERAL_Q(2,:),'b','LineWidth',linW);
plot(0:1:Horizon-1,SVAR.BP_DM_G_SIMPLE_Q(2,:),'r--','LineWidth',linW);
title('Blanchard-Perotti: Spending Multiplier','FontSize',font_num)
set(gca,'YTick',-0.5:0.5:2.5)
set(gca,'XTick',0:4:Horizon-1)
axis([0 40 -0.75 2.5])
box on

subplot(323)
hline(0,'k');
hold on
a = SVAR.PF_DM_T_GENERAL_Q(1,:);
b = SVAR.PF_DM_T_GENERAL_Q(3,:);
[~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
hold on
plot(0:1:Horizon-1,SVAR.PF_DM_T_GENERAL_Q(2,:),'b','LineWidth',linW);
plot(0:1:Horizon-1,SVAR.PF_DM_T_SIMPLE_Q(2,:),'r--','LineWidth',linW);
title('Penalty Function: Tax Multiplier','FontSize',font_num)
set(gca,'YTick',-0.5:0.5:2.5)
set(gca,'XTick',0:4:Horizon-1)
axis([0 40 -0.75 2.5])
box on

subplot(324)
hline(0,'k');
hold on
a = SVAR.PF_DM_G_GENERAL_Q(1,:);
b = SVAR.PF_DM_G_GENERAL_Q(3,:);
[~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
hold on
plot(0:1:Horizon-1,SVAR.PF_DM_G_GENERAL_Q(2,:),'b','LineWidth',linW);
plot(0:1:Horizon-1,SVAR.PF_DM_G_SIMPLE_Q(2,:),'r--','LineWidth',linW);
title('Penalty Function: Spending Multiplier','FontSize',font_num)
set(gca,'YTick',-0.5:0.5:2.5)
set(gca,'XTick',0:4:Horizon-1)
axis([0 40 -0.75 2.5])
box on

subplot(325)
hline(0,'k');
hold on
a = SVAR.PROXY_DM_T_GENERAL_Q(1,:);
b = SVAR.PROXY_DM_T_GENERAL_Q(3,:);
[~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
hold on
plot(0:1:Horizon-1,SVAR.PROXY_DM_T_GENERAL_Q(2,:),'b','LineWidth',linW);
plot(0:1:Horizon-1,SVAR.PROXY_DM_T_SIMPLE_Q(2,:),'r--','LineWidth',linW);
title('Proxy SVAR: Tax Multiplier','FontSize',font_num)
set(gca,'YTick',-0.5:0.5:2.5)
set(gca,'XTick',0:4:Horizon-1)
axis([0 40 -0.75 2.5])
box on

lab1 = sprintf('General Fiscal Rule');
lab2 = sprintf('Simple Fiscal Rule{   }');
     
legend([h1,h2],{lab1,lab2},'Orientation','Horizontal','Location',[0.35,0,0.3,0.05],'FontSize',10,'box','off')

dim = [12,8];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(fig,'-depsc','./results/Figure2');
print(fig,'-dpdf' ,'./results/Figure2');

%% Table 1

str1 = '\\multicolumn{1}{L{8cm}}{$%20.20s$}   & %5.2f             & %5.2f           & %5.2f           \\\\ \n';
str2 = '\\multicolumn{1}{L{8cm}}{$%20.20s$}   & [%5.2f~~~%5.2f]   & [%5.2f~~~%5.2f] & [%5.2f~~~%5.2f] \\\\ \n';
str3 = '\\multicolumn{1}{L{8cm}}{$%20.20s$}   &                   & [%5.2f~~~%5.2f] & [%5.2f~~~%5.2f] \\\\ \n';
str4 = '\\multicolumn{1}{L{8cm}}{$%20.20s$}   & %5.2f             & %5.2f           &                 \\\\ \n';
str5 = '\\multicolumn{1}{L{8cm}}{$%20.20s$}   &                   & [%5.2f~~~%5.2f] &                 \\\\ \n';
str6 = '\\multicolumn{1}{L{8cm}}{$%20.20s$}   &                   &                 &                 \\\\ \n';

%  
fprintf(' \n \\multicolumn{4}{l}{\\rule[-1.5mm]{0mm}{6mm}(A.) \\textit{Tax Rule}}       \\ \\midrule \n')
fprintf(str1, '\psi^{tr}_{0,gdp}',BP_EETA_GENERAL(1,1), SVAR.PF_EETAT_GENERAL_Q(2,3)    ,SVAR.PROXY_EETA_GENERAL_Q(2,2))
fprintf(str3, '                 '                     , SVAR.PF_EETAT_GENERAL_Q([1 3],3), SVAR.PROXY_EETA_GENERAL_Q([1 3],2))
fprintf(str1, '\psi^{tr}_{0,\pi}',BP_EETA_GENERAL(2,1), SVAR.PF_EETAT_GENERAL_Q(2,4)    ,SVAR.PROXY_EETA_GENERAL_Q(2,3))
fprintf(str3, '                 '                     , SVAR.PF_EETAT_GENERAL_Q([1 3],4), SVAR.PROXY_EETA_GENERAL_Q([1 3],3))
fprintf(str1, '\psi^{tr}_{0,r}' ,BP_EETA_GENERAL(3,1) , SVAR.PF_EETAT_GENERAL_Q(2,5)    ,SVAR.PROXY_EETA_GENERAL_Q(2,4))
fprintf(str3, '                 '                     , SVAR.PF_EETAT_GENERAL_Q([1 3],5), SVAR.PROXY_EETA_GENERAL_Q([1 3],4))
fprintf(str1, '\psi^{tr}_{0,g}',SVAR.BP_EETA_TG_Q(2)  , SVAR.PF_EETAT_GENERAL_Q(2,2)    ,SVAR.PROXY_EETA_GENERAL_Q(2,1))
fprintf(str2, '               ',SVAR.BP_EETA_TG_Q([1 3])   , SVAR.PF_EETAT_GENERAL_Q([1 3],2), SVAR.PROXY_EETA_GENERAL_Q([1 3],1))
fprintf(' \n \\multicolumn{4}{l}{\\rule[-1.5mm]{0mm}{6mm}(B.) \\textit{Government Spending Rule}}       \\ \\midrule \n')
fprintf(str4, '\psi^{tr}_{0,gdp}',BP_EETA_GENERAL(1,2), SVAR.PF_EETAG_GENERAL_Q(2,3)    )
fprintf(str5, '                 '                     , SVAR.PF_EETAG_GENERAL_Q([1 3],3))
fprintf(str4, '\psi^{tr}_{0,\pi}',BP_EETA_GENERAL(2,2), SVAR.PF_EETAG_GENERAL_Q(2,4))
fprintf(str5, '                 '                     , SVAR.PF_EETAG_GENERAL_Q([1 3],4))
fprintf(str4, '\psi^{tr}_{0,r}' ,BP_EETA_GENERAL(3,2) , SVAR.PF_EETAG_GENERAL_Q(2,5))
fprintf(str5, '                 '                     , SVAR.PF_EETAG_GENERAL_Q([1 3],5))
fprintf(str4, '\psi^{tr}_{0,tr}',0.00                 , SVAR.PF_EETAG_GENERAL_Q(2,1))
fprintf(str6, '                 '                     )





