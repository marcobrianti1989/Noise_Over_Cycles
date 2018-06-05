tic

%-----------------
% Housekeeping    
%----------------

clear all;  
close all;
clc; 

addpath('./results')
addpath('./auxfiles')

transp90 = 0.25; % Transparency of shaded area in figures
font_num = 12;   % Font size in figures
linW = 2;        % Line width in figures

%--------------------
% Model Specification  
%--------------------
do_Figure = 3;   % Select 3 to 7 to produce the associated figure number
                 % NOTE 1: to produce Figures 4 to 7 you first need to
                 % generate the output for figure 3
                 % NOTE 2: do_Figure = 3 generates Tables 2 and 3
                 

i_var_instr_simple  = {'DTFP_UTIL'};                         % Instrument for proxy SVAR under simple rule
i_var_instr_general = {'DTFP_UTIL','RESID08','HAMILTON3YP'}; % Instrument for proxy SVAR under general rule

switch do_Figure
    case 3
        VAR.select_vars = {'GDP_S','TB3MS','CPI_PIQ4','G_S','TAX_S'};
    case 4
        VAR.select_vars = {'GDP_S','TB3MS','CPI_PIQ4','G_S','TAX_S','MUNI1Y','PDVMILY'};
    case 5
        VAR.select_vars = {'GDP','TB3MS','CPI_PIQ4','G','TAX'};
    case 6
        VAR.select_vars = {'GDP_S','TB3MS','CPI_PIQ4','TAX_S','G_S'};
    case 7
        VAR.select_vars = {'GDP_S','TB3MS','CPI_PIQ4','G_S','TAX_S'};
end

p               = 4;      % Number of lags
nex_            = 1;      % Constant
n               = size(VAR.select_vars,2); % Number of variables
T0              = p;      % Size of pre-sample for Minnesota Prior
Horizon         = 41;     % Horizon for Impulse Responses
nd              = 20000;  % Number of draws in MC chain
bburn           = 0.2*nd; % Share of draws to burn

ptileVEC = [0.16 0.50 0.84]; % Percentiles for posterior distributions

%--------------------
% Load Data  
%--------------------

load DATASET
GDPRATIO = DATASET.GDPRATIO;

YY                   = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,VAR.select_vars)));
VAR.ivar_str         = DATASET.FIGLABELS(cell2mat(values(DATASET.MAP,VAR.select_vars)),1);
VAR.MAP              = containers.Map([VAR.select_vars],1:size(YY,2));
VAR.mm_simple        = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,i_var_instr_simple)));
VAR.mm_general       = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,i_var_instr_general)));
VAR.mF               = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,{'TAXNARRATIVE','PDVMILY'}))); % For F-tests

if do_Figure == 4;
    YY = YY(13:end,:);
    VAR.mm_simple = VAR.mm_simple(13:end,:);
    VAR.mm_general = VAR.mm_general(13:end,:);
end


M_simple             = VAR.mm_simple(T0+1:end,:);
nIV_simple           = size(M_simple,2);

M_general_full       = VAR.mm_general(T0+1:end,[1 3]);
TIV = size(VAR.mm_general,1)-T0;
[row, col] = find(isnan(VAR.mm_general));
M_general_short      = VAR.mm_general(max(row)+1:end,2);

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


%-------------------------
% Run F-tests
%-------------------------
if do_Figure == 3;

    MF                   = [ones(TIV,1) VAR.mF(T0+1:end,:)];
    MF2                  = [ones(TIV-max(row)+T0,1) VAR.mF(max(row)+1:end,:)];
    % Regressions for relevance
    %--------------------------
    % GDP VAR residual
    resultu_gdp = ols(U(ndum+1:end,1),[ones(size(M_general_full,1),1) M_general_full]);
    resultr = ols(U(ndum+1:end,1),ones(size(M_general_full,1),1));
    [fstat_gdp, ~] = waldf(resultr,resultu_gdp);
    % Inflation VAR residual
    resultu_inf = ols(U(ndum+1:end,3),[ones(size(M_general_full,1),1) M_general_full]);
    resultr = ols(U(ndum+1:end,3),ones(size(M_general_full,1),1));
    [fstat_inf, ~] = waldf(resultr,resultu_inf);
    % Interest rate VAR residual
    resultu_rate = ols(U(max(row)+ndum+1-T0:end,2),[ones(size(M_general_short,1),1) M_general_short]);
    resultr = ols(U(max(row)+ndum+1-T0:end,2),ones(size(M_general_short,1),1));
    [fstat_rate, ~] = waldf(resultr,resultu_rate);
    
    % Regressions for exogeneity
    %---------------------------
    
    % TFP utilization adjusted
    resultu_tfp = ols(M_general_full(:,1),MF);
    resultr = ols(M_general_full(:,1),ones(size(M_general_full,1),1));
    [fstat_tfp, ~] = waldf(resultr,resultu_tfp);
    
    % Oil shocks
    resultu_oil = ols(M_general_full(:,2),MF);
    resultr = ols(M_general_full(:,2),ones(size(M_general_full,1),1));
    [fstat_oil, ~] = waldf(resultr,resultu_oil);
    
    % Monetary policy shocks
    resultu_mp = ols(M_general_short(:,1),MF2);
    resultr = ols(M_general_short(:,1),ones(size(M_general_short,1),1));
    [fstat_mp, ~] = waldf(resultr,resultu_mp);
end


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
NF_DM_T_SIMPLE = zeros(nd-bburn,Horizon);
NF_DM_G_SIMPLE = zeros(nd-bburn,Horizon);
NF_DM_T_GENERAL = zeros(nd-bburn,Horizon);
NF_DM_G_GENERAL = zeros(nd-bburn,Horizon);
NF_EETAT_GENERAL = zeros(nd-bburn,n,1);
NF_EETAG_GENERAL = zeros(nd-bburn,n,1);
NF_EETAT_SIMPLE = zeros(nd-bburn,n,1);
NF_EETAG_SIMPLE = zeros(nd-bburn,n,1);

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
        
    Bdraw = reshape(Bdraw,n*p+nex_,n); % Reshape Bdraw from vector to matrix
    Udraw = Y-X*Bdraw;                 % Store residuals for IV regressions
    
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
                       
        [PROXY_EETA, PROXY_IRF_SIMPLE] =  proxy_svar_nonfiscal(Sigmadraw,Udraw(ndum+1:end,:),n,p,F,J,Horizon,M_simple,[],GDPRATIO,do_Figure);
        NF_EETAG_SIMPLE(record-bburn,:) = PROXY_EETA(:,1);
        NF_EETAT_SIMPLE(record-bburn,:) = PROXY_EETA(:,2);
        NF_DM_G_SIMPLE(record-bburn,:)  = PROXY_IRF_SIMPLE(:,1);
        NF_DM_T_SIMPLE(record-bburn,:)  = PROXY_IRF_SIMPLE(:,2);
        
        [PROXY_EETA, PROXY_IRF_GENERAL] =  proxy_svar_nonfiscal(Sigmadraw,Udraw(ndum+1:end,:),n,p,F,J,Horizon,M_general_full,M_general_short,GDPRATIO,do_Figure);
        NF_EETAG_GENERAL(record-bburn,:) = PROXY_EETA(:,1);
        NF_EETAT_GENERAL(record-bburn,:) = PROXY_EETA(:,2);
        NF_DM_G_GENERAL(record-bburn,:)  = PROXY_IRF_GENERAL(:,1);
        NF_DM_T_GENERAL(record-bburn,:)  = PROXY_IRF_GENERAL(:,2);
    end
end

SVAR.NF_DM_T_SIMPLE_Q     = quantile(NF_DM_T_SIMPLE,ptileVEC);
SVAR.NF_DM_G_SIMPLE_Q     = quantile(NF_DM_G_SIMPLE,ptileVEC);
SVAR.NF_DM_T_GENERAL_Q    = quantile(NF_DM_T_GENERAL,ptileVEC);
SVAR.NF_DM_G_GENERAL_Q    = quantile(NF_DM_G_GENERAL,ptileVEC);
SVAR.NF_EETAT_GENERAL_Q   = quantile(NF_EETAT_GENERAL,ptileVEC);
SVAR.NF_EETAG_GENERAL_Q   = quantile(NF_EETAG_GENERAL,ptileVEC);
SVAR.NF_EETAT_SIMPLE_Q   = quantile(NF_EETAT_SIMPLE,ptileVEC);
SVAR.NF_EETAG_SIMPLE_Q   = quantile(NF_EETAG_SIMPLE,ptileVEC);


switch do_Figure
    case 3
        savefileirf = strcat('./results/Results_Figure3','.mat');
    case 4
        savefileirf = strcat('./results/Results_Figure4','.mat');
    case 5
        savefileirf = strcat('./results/Results_Figure5','.mat');
    case 6
        savefileirf = strcat('./results/Results_Figure6','.mat');
    case 7
        savefileirf = strcat('./results/Results_Figure7','.mat');
end

save(savefileirf,'SVAR'); 

toc

%% FIGURE 3

switch do_Figure
    case 3
        fig = figure(3);
        subplot(121)
        hline(0,'k');
        hold on
        a = SVAR.NF_DM_T_GENERAL_Q(1,:);
        b = SVAR.NF_DM_T_GENERAL_Q(3,:);
        [~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
        hold on
        h1 = plot(0:1:Horizon-1,SVAR.NF_DM_T_GENERAL_Q(2,:),'b','LineWidth',linW);
        h2 = plot(0:1:Horizon-1,SVAR.NF_DM_T_SIMPLE_Q(2,:),'r--','LineWidth',linW);
        title('Tax Multiplier','FontSize',font_num)
        set(gca,'YTick',-0.5:0.5:2.5)
        set(gca,'XTick',0:4:Horizon-1)
        axis([0 40 -0.75 2.5])
        box on

        subplot(122)
        hline(0,'k');
        hold on
        a = SVAR.NF_DM_G_GENERAL_Q(1,:);
        b = SVAR.NF_DM_G_GENERAL_Q(3,:);
        [~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
        hold on
        plot(0:1:Horizon-1,SVAR.NF_DM_G_GENERAL_Q (2,:),'b','LineWidth',linW);
        plot(0:1:Horizon-1,SVAR.NF_DM_G_SIMPLE_Q (2,:),'r--','LineWidth',linW);
        title('Spending Multiplier','FontSize',font_num)
        set(gca,'YTick',-0.5:0.5:2.5)
        set(gca,'XTick',0:4:Horizon-1)
        axis([0 40 -0.75 2.5])
        box on

        lab1 = sprintf('General Fiscal Rule{   }');
        lab2 = sprintf('Simple Fiscal Rule');

        legend([h1,h2],{lab1,lab2},'Orientation','Horizontal','Location',[0.35,0,0.3,0.05],'FontSize',10,'box','off')

        dim = [12,4];
        set(gcf,'paperpositionmode','manual','paperunits','inches');
        set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
        print(fig,'-dpdf' ,strcat('./results/Figure3'));
        print(fig,'-depsc',strcat('./results/Figure3'));

        %% TABLE 2

        str1 = '\\multicolumn{1}{L{9cm}}{$%20.20s$}   & %5.2f     & %5.2f   &         \\\\ \n';
        str2 = '\\multicolumn{1}{L{9cm}}{$%20.20s$}   & [%5.2f]   & [%5.2f] &         \\\\ \n';
        str3 = '\\multicolumn{1}{L{9cm}}{$%20.20s$}   &           &         & %5.2f   \\\\ \n';
        str4 = '\\multicolumn{1}{L{9cm}}{$%20.20s$}   &           &         & [%5.2f] \\\\ \n';
        str5 = '\\multicolumn{1}{L{9cm}}{$%20.20s$}   & %5.2f     & %5.2f   & %5.2f   \\\\ \n';
        str6 = '\\multicolumn{1}{L{9cm}}{$%20.20s$}   & [%5.2f]   & [%5.2f] & [%5.2f] \\\\ \n';

        fprintf(' \n \\multicolumn{4}{l}{\\rule[-1.5mm]{0mm}{6mm}(A.) \\textit{Relevance of Non-Fiscal Proxies}} \\ \\midrule \n')
        fprintf(str1, 'm_{tfp}',100*resultu_gdp.beta(2), 100*resultu_inf.beta(2))
        fprintf(str2, '       ',100*resultu_gdp.bstd(2), 100*resultu_inf.bstd(2))
        fprintf(str1, 'm_{oil}',100*resultu_gdp.beta(3), 100*resultu_inf.beta(3))
        fprintf(str2, '       ',100*resultu_gdp.bstd(3), 100*resultu_inf.bstd(3))
        fprintf(str3, 'm_{r}'  ,100*resultu_rate.beta(2)                        )
        fprintf(str4, '       ',100*resultu_rate.bstd(2)                        )
        fprintf(str5, '       ',fstat_gdp, fstat_inf, fstat_rate                )
        fprintf(' \n \\multicolumn{4}{l}{\\rule[-1.5mm]{0mm}{6mm}(A.) \\textit{Exogeneity of Non-Fiscal Proxies}} \\ \\midrule \n')
        fprintf(str5, 'm_{tax}',resultu_tfp.beta(2), resultu_oil.beta(2), resultu_mp.beta(2))
        fprintf(str6, '       ',resultu_tfp.bstd(2), resultu_oil.bstd(2), resultu_mp.bstd(2))
        fprintf(str5, 'm_{g}'  ,resultu_tfp.beta(3)/100, resultu_oil.beta(3)/100, resultu_mp.beta(3)/100)
        fprintf(str6, '       ',resultu_tfp.bstd(3)/100, resultu_oil.bstd(3)/100, resultu_mp.bstd(3)/100)
        fprintf(str5, '       ',fstat_tfp, fstat_oil, fstat_mp                )

        
        %% TABLE 3

        str1 = '\\multicolumn{1}{L{9cm}}{$%20.20s$}   & %5.2f             & %5.2f            \\\\ \n';
        str2 = '\\multicolumn{1}{L{9cm}}{$%20.20s$}   & [%5.2f~~~%5.2f]   & [%5.2f~~~%5.2f]  \\\\ \n';
        str3 = '\\multicolumn{1}{L{9cm}}{$%20.20s$}   & %5.2f             &                  \\\\ \n';
        str4 = '\\multicolumn{1}{L{9cm}}{$%20.20s$}   & [%5.2f~~~%5.2f]   &                  \\\\ \n';
        str5 = '\\multicolumn{1}{L{9cm}}{$%20.20s$}   &                   &                  \\\\ \n';

        fprintf(' \n \\multicolumn{4}{l}{\\rule[-1.5mm]{0mm}{6mm}(A.) \\textit{Tax Rule}}       \\ \\midrule \n')
        fprintf(str1, '\psi^{tr}_{0,gdp}',SVAR.NF_EETAT_GENERAL_Q(2,1)    , SVAR.NF_EETAT_SIMPLE_Q(2,1)    )
        fprintf(str2, '                 ',SVAR.NF_EETAT_GENERAL_Q([1 3],1), SVAR.NF_EETAT_SIMPLE_Q([1 3],1))
        fprintf(str3, '\psi^{tr}_{0,\pi}',SVAR.NF_EETAT_GENERAL_Q(2,3)    )
        fprintf(str4, '                 ',SVAR.NF_EETAT_GENERAL_Q([1 3],3))
        fprintf(str3, '\psi^{tr}_{0,r}'  ,SVAR.NF_EETAT_GENERAL_Q(2,2)    )
        fprintf(str4, '                 ',SVAR.NF_EETAT_GENERAL_Q([1 3],2))
        fprintf(str3, '\psi^{tr}_{0,g}'  ,SVAR.NF_EETAT_GENERAL_Q(2,4)    )
        fprintf(str4, '                 ',SVAR.NF_EETAT_GENERAL_Q([1 3],4))
        fprintf(' \n \\multicolumn{4}{l}{\\rule[-1.5mm]{0mm}{6mm}(B.) \\textit{Government Spending Rule}}       \\ \\midrule \n')
        fprintf(str1, '\psi^{g}_{0,gdp}' ,SVAR.NF_EETAG_GENERAL_Q(2,1)    , SVAR.NF_EETAG_SIMPLE_Q(2,1)    )
        fprintf(str2, '                 ',SVAR.NF_EETAG_GENERAL_Q([1 3],1), SVAR.NF_EETAG_SIMPLE_Q([1 3],1))
        fprintf(str3, '\psi^{g}_{0,\pi}' ,SVAR.NF_EETAG_GENERAL_Q(2,3)    )
        fprintf(str4, '                 ',SVAR.NF_EETAG_GENERAL_Q([1 3],3))
        fprintf(str3, '\psi^{g}_{0,r}'   ,SVAR.NF_EETAG_GENERAL_Q(2,2)    )
        fprintf(str4, '                 ',SVAR.NF_EETAG_GENERAL_Q([1 3],2))
        fprintf(str3, '\psi^{g}_{0,tr}'  ,SVAR.NF_EETAG_GENERAL_Q(2,5)    )
        fprintf(str5, '                 ')
    case 4
        BASE = load('./results/Results_Figure3');
        fig = figure(4);
        subplot(121)
        hline(0,'k');
        hold on
        a = SVAR.NF_DM_T_GENERAL_Q(1,:);
        b = SVAR.NF_DM_T_GENERAL_Q(3,:);
        [~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
        hold on
        h1 = plot(0:1:Horizon-1,SVAR.NF_DM_T_GENERAL_Q(2,:),'b','LineWidth',linW);
        h2 = plot(0:1:Horizon-1,BASE.SVAR.NF_DM_T_GENERAL_Q(2,:),'r--','LineWidth',linW);
        title('Tax Multiplier','FontSize',font_num)
        set(gca,'YTick',-0.5:0.5:2.5)
        set(gca,'XTick',0:4:Horizon-1)
        axis([0 40 -0.75 2.5])
        box on

        subplot(122)
        hline(0,'k');
        hold on
        a = SVAR.NF_DM_G_GENERAL_Q(1,:);
        b = SVAR.NF_DM_G_GENERAL_Q(3,:);
        [~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
        hold on
        plot(0:1:Horizon-1,SVAR.NF_DM_G_GENERAL_Q (2,:),'b','LineWidth',linW);
        plot(0:1:Horizon-1,BASE.SVAR.NF_DM_G_GENERAL_Q (2,:),'r--','LineWidth',linW);
        title('Spending Multiplier','FontSize',font_num)
        set(gca,'YTick',-0.5:0.5:2.5)
        set(gca,'XTick',0:4:Horizon-1)
        axis([0 40 -0.75 2.5])
        box on

        lab1 = sprintf('Model with Fiscal News{   }');
        lab2 = sprintf('Baseline');

        legend([h1,h2],{lab1,lab2},'Orientation','Horizontal','Location',[0.35,0,0.3,0.05],'FontSize',10,'box','off')

        dim = [12,4];
        set(gcf,'paperpositionmode','manual','paperunits','inches');
        set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
        print(fig,'-dpdf' ,strcat('./results/Figure4'));
        print(fig,'-depsc',strcat('./results/Figure4'));
    case 5
        BASE = load('./results/Results_Figure3');
        fig = figure(5);
        subplot(121)
        hline(0,'k');
        hold on
        a = SVAR.NF_DM_T_GENERAL_Q(1,:);
        b = SVAR.NF_DM_T_GENERAL_Q(3,:);
        [~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
        hold on
        h1 = plot(0:1:Horizon-1,SVAR.NF_DM_T_GENERAL_Q(2,:),'b','LineWidth',linW);
        h2 = plot(0:1:Horizon-1,BASE.SVAR.NF_DM_T_GENERAL_Q(2,:),'r--','LineWidth',linW);
        title('Tax Multiplier','FontSize',font_num)
        set(gca,'YTick',-0.5:0.5:2.5)
        set(gca,'XTick',0:4:Horizon-1)
        axis([0 40 -0.75 2.5])
        box on

        subplot(122)
        hline(0,'k');
        hold on
        a = SVAR.NF_DM_G_GENERAL_Q(1,:);
        b = SVAR.NF_DM_G_GENERAL_Q(3,:);
        [~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
        hold on
        plot(0:1:Horizon-1,SVAR.NF_DM_G_GENERAL_Q (2,:),'b','LineWidth',linW);
        plot(0:1:Horizon-1,BASE.SVAR.NF_DM_G_GENERAL_Q (2,:),'r--','LineWidth',linW);
        title('Spending Multiplier','FontSize',font_num)
        set(gca,'YTick',-0.5:0.5:2.5)
        set(gca,'XTick',0:4:Horizon-1)
        axis([0 40 -0.75 2.5])
        box on

        lab1 = sprintf('Data Not Detrended{   }');
        lab2 = sprintf('Baseline');

        legend([h1,h2],{lab1,lab2},'Orientation','Horizontal','Location',[0.35,0,0.3,0.05],'FontSize',10,'box','off')

        dim = [12,4];
        set(gcf,'paperpositionmode','manual','paperunits','inches');
        set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
        print(fig,'-dpdf' ,strcat('./results/Figure5'));
        print(fig,'-depsc',strcat('./results/Figure5'));
    case 6
        BASE = load('./results/Results_Figure3');
        fig = figure(6);
        subplot(121)
        hline(0,'k');
        hold on
        a = SVAR.NF_DM_T_GENERAL_Q(1,:);
        b = SVAR.NF_DM_T_GENERAL_Q(3,:);
        [~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
        hold on
        h1 = plot(0:1:Horizon-1,SVAR.NF_DM_T_GENERAL_Q(2,:),'b','LineWidth',linW);
        h2 = plot(0:1:Horizon-1,BASE.SVAR.NF_DM_T_GENERAL_Q(2,:),'r--','LineWidth',linW);
        title('Tax Multiplier','FontSize',font_num)
        set(gca,'YTick',-0.5:0.5:2.5)
        set(gca,'XTick',0:4:Horizon-1)
        axis([0 40 -0.75 2.5])
        box on

        subplot(122)
        hline(0,'k');
        hold on
        a = SVAR.NF_DM_G_GENERAL_Q(1,:);
        b = SVAR.NF_DM_G_GENERAL_Q(3,:);
        [~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
        hold on
        plot(0:1:Horizon-1,SVAR.NF_DM_G_GENERAL_Q (2,:),'b','LineWidth',linW);
        plot(0:1:Horizon-1,BASE.SVAR.NF_DM_G_GENERAL_Q (2,:),'r--','LineWidth',linW);
        title('Spending Multiplier','FontSize',font_num)
        set(gca,'YTick',-0.5:0.5:2.5)
        set(gca,'XTick',0:4:Horizon-1)
        axis([0 40 -0.75 2.5])
        box on

        lab1 = sprintf('Gov. Spending Ordered Second{   }');
        lab2 = sprintf('Baseline');

        legend([h1,h2],{lab1,lab2},'Orientation','Horizontal','Location',[0.35,0,0.3,0.05],'FontSize',10,'box','off')

        dim = [12,4];
        set(gcf,'paperpositionmode','manual','paperunits','inches');
        set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
        print(fig,'-dpdf' ,strcat('./results/Figure6'));
        print(fig,'-depsc',strcat('./results/Figure6'));
    case 7
        BASE = load('./results/Results_Figure3');
        fig = figure(7);
        subplot(121)
        hline(0,'k');
        hold on
        a = SVAR.NF_DM_T_GENERAL_Q(1,:);
        b = SVAR.NF_DM_T_GENERAL_Q(3,:);
        [~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
        hold on
        h1 = plot(0:1:Horizon-1,SVAR.NF_DM_T_GENERAL_Q(2,:),'b','LineWidth',linW);
        h2 = plot(0:1:Horizon-1,BASE.SVAR.NF_DM_T_GENERAL_Q(2,:),'r--','LineWidth',linW);
        title('Tax Multiplier','FontSize',font_num)
        set(gca,'YTick',-0.5:0.5:2.5)
        set(gca,'XTick',0:4:Horizon-1)
        axis([0 40 -0.75 2.5])
        box on

        subplot(122)
        hline(0,'k');
        hold on
        a = SVAR.NF_DM_G_GENERAL_Q(1,:);
        b = SVAR.NF_DM_G_GENERAL_Q(3,:);
        [~,~]=jbfill(0:1:Horizon-1,a,b,'b','none',1,transp90);
        hold on
        plot(0:1:Horizon-1,SVAR.NF_DM_G_GENERAL_Q (2,:),'b','LineWidth',linW);
        plot(0:1:Horizon-1,BASE.SVAR.NF_DM_G_GENERAL_Q (2,:),'r--','LineWidth',linW);
        title('Spending Multiplier','FontSize',font_num)
        set(gca,'YTick',-0.5:0.5:2.5)
        set(gca,'XTick',0:4:Horizon-1)
        axis([0 40 -0.75 2.5])
        box on

        lab1 = sprintf('Alternative Scaling of Shocks{   }');
        lab2 = sprintf('Baseline');

        legend([h1,h2],{lab1,lab2},'Orientation','Horizontal','Location',[0.35,0,0.3,0.05],'FontSize',10,'box','off')

        dim = [12,4];
        set(gcf,'paperpositionmode','manual','paperunits','inches');
        set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
        print(fig,'-dpdf' ,strcat('./results/Figure7'));
        print(fig,'-depsc',strcat('./results/Figure7'));
end





