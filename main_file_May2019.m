%*************************************************************************%
%    Main File, last change 04/27/2019
%    Code by Brianti, Marco e Cormun, Vito
%*************************************************************************%

clear
close all

create_TFP_surprise;
clear
load('TechShock_identification.mat');

% Technical Parameters
lags                = 4;             % Number of lags in the first step (deriving Ztilde)
leads               = 0;             % Number of leads in the first step (deriving Ztilde)
H                   = 20;            % IRFs horizon
lags_LP             = 2;            % Number of lags in the Local Projection
which_trend         = 'quad';        %'BP', 'HP', 'lin', 'quad', 'none', 'demean' for Local Projection
which_Z             = {'1','2','3','4','5'}; % Which Forecast Revision: RGDP, NGDP, RCONS, INDPROD, RINV. If it is more than one it takes the first PC
which_shock         = {'Sentiment','Tech'};      % Tech, News, Sentiment
loc_start_exogenous = 0;             % Exogenous start
diff_LP             = 0;             % LP in levels or differences
nPC_first           = 3;             % Number of Principal Components in the first stage
nPC_LP              = 2;             % Number of Principal Components in the second stage
norm_SHOCK          = 0;             % Divide shock over its own variance
printIRFs           = 1;             % Print IRFs
printVD             = 0;              % Print Variance Decompositions
nsimul              = 500;           % number of simulations for bootstrap
control_pop         = 0;             % Divide GDP, Cons, Hours, Investment over population
varlist             = {'ProfitAdj','ProfitaftTax','Dividends','RealGDP','Cashflow'};%{'RealGDP','RealInvestment','RealCons','HoursAll'};%{'ProfitAdj','ProfitbefTax','ProfitaftTax','Dividends','RealGDP','Cashflow'};%{'HoursAll','Hours','HoursPerPerson','UnempRate','Employment'};%,'GDPDefl','FFR'}; % Define endogenous variables for LP
varlist_graph       = varlist; %{'RealGDP','RealInvestment','RealCons','Hours'};

% Read main dataset
filename                    = 'main_file';
sheet                       = 'Sheet1';
range                       = 'B1:EG300';
do_truncation               = 0; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data2(filename, sheet, range, do_truncation);
dataset                     = [dataset; NaN(leads,size(dataset,2))]; % Adding some NaN at the end for technical reasons

% Assess name to each variable
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

%*************************************************************************%
%                                                                         %
%                     1st stage - Deriving sentiment                      %
%                                                                         %
%*************************************************************************%

%Building Zt - Forecast Revisions from SPF and Michigan Index
create_Z;

% Define Variables
[TFP_trunc, trunc1, trunc2] = truncate_data(TFP);
TFPBP                       = bpass(TFP_trunc,4,32);
TFPBP                       = [TFPBP; NaN(length(TFP) - length(TFPBP),1)];
dTFP                        = [NaN; diff(TFP)];
PC                          = [PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9];
PC_first                    = PC(:,1:nPC_first);
SHOCKS_NARRATIVE            = [ TAXNARRATIVE  PDVMILY HAMILTON3YP RESID08 ];  % RESID08 reduce number of shocks

% Loop over different surveys
for iw = 1:length(which_Z)
      % Define dependent variable
      eval(['Z = Z', which_Z{iw},';']);
      Yiw(:,iw)                      = Z;
end
% Get first PC
[Y_cut, tt, tt2]   = truncate_data(Yiw);
YPCs               = get_principal_components(Y_cut);
YPCs               = [NaN(tt-1,size(Yiw,2)); YPCs; NaN(size(Yiw,1)-tt2,size(Yiw,2))];
Y                  = YPCs(:,1);
% Define Regressors and Dependent Variable
X_contemporaneous           = [TFPBP];% SHOCKS_NARRATIVE];% [TFPBP SHOCKS_NARRATIVE]; %[TFPBP];%
X_lag                       = [TFPBP PC_first]; % SHOCKS_NARRATIVE];%[TFPBP PC SHOCKS_NARRATIVE]; %[TFPBP PC];%
X_lead                      = TFPBP;
% Control Regression
[~, Zhat, Ztilde, regressor] = lead_lag_matrix_regression(Y,X_lead,...
      leads,X_lag,lags,X_contemporaneous);
%Ztilde                       = (nanvar(Ztilde)^0.5)*randn(length(Ztilde),1);

%*************************************************************************%
%                                                                         %
%                     2nd stage - Local Projections                       %
%                                                                         %
%*************************************************************************%

% Create Inflation from Price Indexes
ninfl           = 4;
PCEInflation    = create_inflation(PriceCPE,ninfl);
CPIInflation    = create_inflation(CPIInflation,ninfl);
CPIDurables     = create_inflation(CPIDurables,ninfl);
CPINonDurables  = create_inflation(CPINonDurables,ninfl);
CPIServices     = create_inflation(CPIServices,ninfl);

% Per capita adjustment
if control_pop == 1
      RealGDP                 = RealGDP - Population;
      RealCons                = RServiceCons + RNonDurableCons - Population;
      RealInvestment          = RealInvestment - Population;
      HoursPerPerson          = HoursAll - Population;
      RealInventories         = RealInventories - Population;
      RealSales               = RealSales - Population;
else
      HoursPerPerson          = HoursAll - Population;
      RealCons                = RServiceCons + RNonDurableCons;
      RealInvestment          = RealInvestment; % + RDurableCons;
end

% Other Transformations
SP500                   = SP500 - GDPDefl;
SpreadBond              = MoodySpreadBaa - MoodySpreadAaa;
TermYield               = TenYTreasury - ThreeMTreasury;
%SpreadBankRate          = LoanPrime - TenYTreasury;
BankLeverage            = BanksTotLiabilities - BanksTotAssets;
CorpEqui2Assets         = NonFinEquity - NonFinTotAssets;
CorpDebt2Assets         = NonFinDebtSecurities - NonFinTotAssets;
ProfitAdj               = CorpProfitsAdj - GDPDefl;
ProfitbefTax            = CorpProfitsbefTax - GDPDefl;
ProfitaftTax            = CorpProfitsNoAdj - GDPDefl;
Dividends               = Dividends - GDPDefl;
UndistProf              = UndistributedProfits - GDPDefl;
Cashflow                = CashFlow - GDPDefl;
CorpEqui2Assets         = NonFinEquity - GDPDefl;
Epay                    = 100*EPayout./exp(RealGDP);%./exp(ValueAddedNFCorp); %./exp(RealGDP);
Drep                    = 100*DebtRep./exp(RealGDP);%./exp(ValueAddedNFCorp); %./exp(RealGDP); 
% ProfAdj                 = CorpProfitsAdj - NonFinEquity;
% ProfbT                  = CorpProfitsbefTax - NonFinEquity;
% Prof                    = CorpProfitsNoAdj - NonFinEquity;
% Div                     = Dividends - NonFinEquity;
% UndistProf              = UndistributedProfits - NonFinEquity;
% Cashflow                = CashFlow - NonFinEquity;
% CorpEqui2Assets         = NonFinEquity - NonFinEquity;
% CorpDebt2Assets         = NonFinDebtSecurities - NonFinTotAssets;
% ProfAdj                 = CorpProfitsAdj - NonFinTotAssets;
% ProfbT                  = CorpProfitsbefTax - NonFinTotAssets;
% Prof                    = CorpProfitsNoAdj - NonFinTotAssets;
% Div                     = Dividends - NonFinTotAssets;
% UndistProf              = UndistributedProfits - NonFinTotAssets;
% Cashflow                = CashFlow - NonFinTotAssets;

% corr(Epay(1:end),RealGDP,'rows','complete')
% corr(Drep(1:end),RealGDP,'rows','complete')



% Matrix of dependent variables - All the variables are in log levels
for i = 1:length(varlist)
      dep_var(:,i) = eval(varlist{i});
end

% Choose the type of transformation
if diff_LP == 1
      dep_var = [nan(1,size(dep_var,2)); diff(dep_var)];
end



% Align timing with SHOCK
dep_var          = dep_var(1+lags:end-leads,:);
PC_LP            = PC(1+lags:end-leads,1:nPC_LP);
Time             = Time(1+lags:end-leads);

% Which shock to plot
for is = 1:length(which_shock)
      if strcmp(which_shock{is},'Sentiment') == 1
            SHOCK = Ztilde;%(1+lags:end-leads); %XXX: BE CAREFUL ZTILDE IS FROM main_file.xls (PREDETERMINED)
            fprintf('\n')
            disp('Sentiment Shock')
            fprintf('\n')
      elseif strcmp(which_shock{is},'Tech') == 1
            SHOCK = TechSurprise(1+lags:end-leads);      %UnantTFPshock(1+lags:end-leads);
            fprintf('\n')
            disp('Unexpected Technology Shock')
            fprintf('\n')
      elseif strcmp(which_shock{is},'News') == 1
            SHOCK = BarskySimsNews(1+lags:end-leads);
            fprintf('\n')
            disp('Barsky and Sims News Shock')
            fprintf('\n')
      end
      IRFname = ['OO_IRFs_',char(which_shock(is)),'Shock_lags',num2str(lags),'_','_leads',num2str(leads),'_lagsLP',num2str(lags_LP),'_trend',which_trend,'_diff',num2str(diff_LP),'_nPC',num2str(nPC_LP),'.png'];
      VDname  = ['VD_',char(which_shock(is)),'Shock_lags',num2str(lags),'_','_leads',num2str(leads),'_lagsLP',num2str(lags_LP),'_trend',which_trend,'_diff',num2str(diff_LP),'_nPC',num2str(nPC_LP)','.png'];
      
      % Normilize Variance of SHOCK
      %       if norm_SHOCK == 1
      %             sdSHOCK          = nanstd(SHOCK);
      %             SHOCK            = SHOCK/sdSHOCK;
      %       end
      
      %Initializating the loop
      for kk = 1:size(dep_var,2)
            % Define inputs for local_projection
            varnamekk                   = varlist{kk};
            disp(['  - Projecting ',varnamekk])
            fprintf('\n')
            depvarkk                           = dep_var(:,kk);
            [~, loc_start(kk,is), loc_end]     = truncate_data([depvarkk SHOCK PC_LP]);
            depvarkk                           = depvarkk(loc_start(kk,is):loc_end);
            SHOCKkk                            = SHOCK(loc_start(kk,is):loc_end);
            pckk                               = PC_LP(loc_start(kk,is):loc_end,:);
            
            % Run local_projection
            [IR{kk},res{kk},tuple{kk},VD{kk},DF{kk},nREG{kk}] = ...
                  local_projection(depvarkk,pckk,SHOCKkk,lags_LP,H,which_trend);
            if diff_LP == 0
                  IRF(kk,:,is) = IR{kk};
            else
                  IRF(kk,:,is) = cumsum(IR{kk});
            end
            VDkk(kk,:)        = VD{kk};
            DFkk(kk,:,is)     = DF{kk};
            nREGkk(kk,:,is)   = nREG{kk};
            % Initiate bootstrap
            tuplekk        = tuple{kk};
            for hh = 1:H
                  tuplekkhh = tuplekk{hh}; % Fix a specific horizon
                  Y                                = tuplekkhh(:,1);
                  X                                = tuplekkhh(:,2:end);
                  [Yboot, Xboot]                   = bb_bootstrap_LP(Y,X,nsimul,lags_LP);
                  for isimul = 1:nsimul
                        B                          = Xboot(:,:,isimul)'*Xboot(:,:,isimul)\...
                              (Xboot(:,:,isimul)'*Yboot(:,isimul));
                        IRF_boot(kk,hh,isimul,is)  = B(1);
                  end
            end
      end
      
      % Select upper and lower bands
      for kk = 1:size(dep_var,2)
            IRF_bootkk = IRF_boot(kk,:,:);
            if diff_LP == 0
                  IRF_boot(kk,:,:)  = IRF_bootkk;
            else
                  IRF_boot(kk,:,:)  = cumsum(IRF_bootkk,2);
            end
      end
      sig              = 0.05;
      sig2             = 0.16;
      for j = 1:size(dep_var,2)
            IRF_up(j,:)   = quantile(squeeze(IRF_boot(j,:,:,is))',1-sig);
            IRF_up2(j,:)  = quantile(squeeze(IRF_boot(j,:,:,is))',1-sig2);
            IRF_low(j,:)  = quantile(squeeze(IRF_boot(j,:,:,is))',sig);
            IRF_low2(j,:) = quantile(squeeze(IRF_boot(j,:,:,is))',sig2);
      end
      
      % Plot IRFs
      plot_IRF(varlist_graph,IRF_low,IRF_low2,IRF_up,IRF_up2,IRF(:,:,is),H,printIRFs,IRFname); %change this function
      % Plot VD
      plot_IRF(varlist,VDkk,VDkk,VDkk,VDkk,VDkk,H,printVD,VDname)
      
end

% Technical Parameters
tech_info_table;

asd

%*************************************************************************%
%                                                                         %
%                        Test Cyclicality - Canova                        %
%                                                                         %
%*************************************************************************%

% Simulate AR(1), Estimate IRFs via LP, Generate Spectral Density
rho             = 0;    % rho = 0, Hnull: flat spectral density, otherwise rho should be estimated from data as im BG
sigm            = 1;    % variance of disturbances (shocks)
T               = 1000; % Asyntotic (Maybe should be the same length) !!!
show_fig_AR_LP  = 1;
show_fig_AR_SD  = 1;
[IRF_AR,IRF_boot_AR] = generate_AR1_LP_IRF(rho,sigm,T,H,nsimul,sig,sig2,show_fig_AR_LP);
% Generate AR(1) Spectral Density
[SD_AR(:,1),SD_AR_UP(:,1),SD_AR_LOW(:,1), ...
      SD_AR_UP2(:,1),SD_AR_LOW2(:,1),SD_AR_MED(:,1),SD_AR_boot,period] ...
      = spectral_density_MA(IRF_AR,IRF_boot_AR,sig,sig2,show_fig_AR_SD);

%Generate Empirical Spectral Density
show_fig_SD  = 1;
for is = 1:length(which_shock)
      disp([which_shock{is}, ' Shock'])
      fprintf('\n')
      for iv = 1:length(varlist)
            disp(['Variable ',varlist{iv}])
            fprintf('\n')
            [SD(:,iv,is),SD_UP(:,iv,is),SD_LOW(:,iv,is), ...
                  SD_UP2(:,iv,is),SD_LOW2(:,iv,is), ...
                  SD_MED(:,iv,is), SD_boot(:,:,iv,is), ~] ...
                  = spectral_density_MA(IRF(iv,:,is),IRF_boot(iv,:,:,is),...
                  sig,sig2,show_fig_SD);
      end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot figure if showfig = 1
per_start = 10;
per_end   = 200;
%plot spectral density of AR(1) from LP
hfig        = findobj('type','figure');
nfig        = length(hfig);
figure(1+nfig)
set(gcf,'Position',[-1919 41 1920 963])
set(gcf,'color','w');
for ii = 1:2
      subplot(1,2,ii)
      plot(period(per_start:per_end)',SD_MED(per_start:per_end,1,ii),'-r','LineWidth',2); hold on; %step dependent
      plot(period(per_start:per_end)',SD_UP2(per_start:per_end,1,ii),'--k','LineWidth',2); hold on;
      plot(period(per_start:per_end)',SD_LOW2(per_start:per_end,1,ii),'--k','LineWidth',2); hold on; %the point estimate is not included in the CI, is it because we don't correct for the bias in the LP?
      plot(period(per_start:per_end)',SD_UP(per_start:per_end,1,ii),'--k','LineWidth',1); hold on;
      plot(period(per_start:per_end)',SD_LOW(per_start:per_end,1,ii),'--k','LineWidth',1); hold on; %the point estimate is not included in the CI, is it because we don't correct for the bias in the LP?
      set(gca,'FontSize',22);
      if ii == 1
            title('Sentiment Shock','fontsize',40)
            xlabel('Periodicity','interpreter','latex','fontsize',30);
            ylabel('Spectral Density','interpreter','latex','fontsize',30);
      elseif ii == 2
            title('Technology Shock','fontsize',40)
      end
      set(gca,'YTickLabel',[]);
      axis tight
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute average spectral density, D1, around the peak  and average
%spectral density around the trough, D2
lpeak_low    = 35; %should be adjusted with steps and IRF horizon
lpeak_up     = 45;
ltrough_low  = 60;
ltrough_up   = 70;
for iv = 1:length(varlist)
      for is = 1:length(which_shock)
            pval(is,iv) = test_canova_peak_sdensity(lpeak_low,lpeak_up,ltrough_low,ltrough_up,...
                  SD_boot(:,:,iv,is),SD_AR,period,nsimul);
      end
end
tech_info_test;
toc

% %Check Overreaction and Underreaction conditional on a shock  - adjust FE
% timing
% if which_Z == '1'
%     [B,Bint] = regress(FE_RGDP(lags+1:end),[Ztilde, ones(length(Ztilde),1)])
% elseif which_Z == '2'
%     [B,Bint] = regress(FE_NGDP(lags+1:end),[Ztilde, ones(length(Ztilde),1)])
% elseif which_Z == '3'
%     [B,Bint] = regress(FE_RC(lags+1:end),[Ztilde, ones(length(Ztilde),1)])
% end
%
% [B,Bint] = regress(FE_RGDP(lags+1:end),[UnantTFPshock(lags+1:end), ones(length(Ztilde),1)])
% [B,Bint] = regress(FE_NGDP(lags+1:end),[UnantTFPshock(lags+1:end), ones(length(Ztilde),1)])
% [B,Bint] = regress(FE_RC(lags+1:end),[UnantTFPshock(lags+1:end), ones(length(Ztilde),1)])
