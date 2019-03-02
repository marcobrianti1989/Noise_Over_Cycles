%*************************************************************************%
%    Main
%
%
%
%    last change 02/14/2019
%
%   todo: Allow for setting up the time period
%   considered, the NFCI, for example, is extremely volatile before 1980
%   Code by Brianti, Marco e Cormun, Vito
%*************************************************************************%
clc
clear
close all
tic
% Start date and end date

% Decide whether or not include policy shocks in the controls

%
% Technical Parameters
lags                = 4;             % Number of lags in the first step (deriving Ztilde)
leads               = 0;             % Number of leads in the first step (deriving Ztilde)
H                   = 20;            % IRFs horizon
lags_LP             = 2;             % Number of lags in the Local Projection
which_trend         = 'quadratic' ;  % BPfilter, HPfilter, linear, quadratic for Local Projection
which_Z             = {'1','2','3','4','5'}; % Which Forecast Revision: RGDP, NGDP, RCONS, INDPROD, RINV. If it is more than one it takes the first PC
which_shock         = {'Sentiment'}; % Tech, News, Sentiment
loc_start_exogenous = 0;       % Exogenous start
diff_LP             = 0;             % LP in levels or differences
nPC_first           = 3;             % Number of Principal Components in the first stage
nPC_LP              = 3;             % Number of Principal Components in the second stage
norm_SHOCK          = 1;             % Divide shock over its own variance
printIRFs           = 1;             % Print IRFs
printVD             = 0;             % Print Variance Decompositions
nsimul              = 2000;           % number of simulations for bootstrap
control_pop         = 0;             % Divide GDP, Cons, Hours, Investment over population
varlist             = {'TermYield','RealGDP','SpreadBond'};%','RealGDP','RealInvestment','RealCons','HoursAll','RealInventories'}; % Define endogenous variables for LP
% 'SpreadBond'  'Leverage'        'ChicagoFedIndex'  'RealExchRate' 'FFR'
% 'SpreadBonds' 'MoodySpreadBaa'  'TermYield'        'FFR'      'Y10Treasury'     'M3Treasury'
% 'RealGDP'     'RealInvestment'  'SpreadBond'       'Leverage' 'ChicagoFedIndex' 'Vix'

% Read main dataset
filename                    = 'main_file';
sheet                       = 'Sheet1';
range                       = 'B1:DW300';
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
SHOCKS_NARRATIVE            = [RESID08]; % TAXNARRATIVE MUNI1Y PDVMILY HAMILTON3YP  reduce number of shocks

% Loop over different surveys
for iw = 1:length(which_Z)
      % Define dependent variable
      eval(['Z = Z', which_Z{iw},';']);
      Y                           = Z;
      % Define Regressors and Dependent Variable
      X_contemporaneous           = [TFPBP];% SHOCKS_NARRATIVE];% [TFPBP SHOCKS_NARRATIVE]; %[TFPBP];%
      X_lag                       = [TFPBP PC_first];% SHOCKS_NARRATIVE];%[TFPBP PC SHOCKS_NARRATIVE]; %[TFPBP PC];%
      X_lead                      = TFPBP;
      % Control Regression
      [~, Zhat, Ztildeiw(:,iw), regressor] = lead_lag_matrix_regression(Y,X_lead,...
            leads,X_lag,lags,X_contemporaneous);
end

% Get first PC over surveys
[Ztilde_cut, tt, tt2]   = truncate_data(Ztildeiw);
ZtildePCs               = get_principal_components(Ztilde_cut);
ZtildePCs               = [NaN(tt-1,size(Ztildeiw,2)); ZtildePCs; NaN(size(Ztildeiw,1)-tt2,size(Ztildeiw,2))];
Ztilde                  = ZtildePCs(:,1);

%*************************************************************************%
%                                                                         %
%                      2nd stage - Local Projections                      %
%                                                                         %
%*************************************************************************%

% Create Inflation from Price Indexes
ninfl = 4;
PCEInflation    = create_inflation(PriceCPE,ninfl);
CPIInflation    = create_inflation(CPIInflation,ninfl);
CPIDurables     = create_inflation(CPIDurables,ninfl);
CPINonDurables  = create_inflation(CPINonDurables,ninfl);
CPIServices     = create_inflation(CPIServices,ninfl);

% Per capita adjustment
if control_pop == 1
      RealGDP                 = RealGDP - Population;
      RealCons                = RealCons - Population;
      RealInvestment          = RealInvestment - Population;
      HoursPerPerson          = HoursAll - Population;
      RealInventories         = RealInventories - Population;
      RealSales               = RealSales - Population;
else
      HoursPerPerson          = HoursAll - Population;
end

% Other Transformations
SP500                   = SP500 - GDPDefl;
SpreadBond              = MoodySpreadBaa - MoodySpreadAaa;
TermYield               = TenYTreasury - ThreeMTreasury;
%SpreadBankRate          = LoanPrime - TenYTreasury;
BankLeverage            = BanksTotLiabilities - BanksTotAssets;
CorpEqui2Assets         = NonFinEquity - NonFinTotAssets;
CorpDebt2Assets         = NonFinDebtSecurities - NonFinTotAssets;

% Matrix of dependen variables - All the variables are in log levels
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
            SHOCK = Ztilde;
            fprintf('\n')
            disp('Sentiment Shock')
            fprintf('\n')
      elseif strcmp(which_shock{is},'Tech') == 1
            SHOCK = UnantTFPshock(1+lags:end-leads);
            fprintf('\n')
            disp('Unexpected Technology Shock')
            fprintf('\n')
      elseif strcmp(which_shock{is},'News') == 1
            SHOCK = BarskySimsNews(1+lags:end-leads);
            fprintf('\n')
            disp('Barsky and Sims News Shock')
            fprintf('\n')
      end
      IRFname = ['IRFs_',char(which_shock(is)),'Shock_lags',num2str(lags),'_','_leads',num2str(leads),'_lagsLP',num2str(lags_LP),'_trend',which_trend,'_diff',num2str(diff_LP),'_nPC',num2str(nPC_LP),'.pdf'];
      VDname  = ['VD_',char(which_shock(is)),'Shock_lags',num2str(lags),'_','_leads',num2str(leads),'_lagsLP',num2str(lags_LP),'_trend',which_trend,'_diff',num2str(diff_LP),'_nPC',num2str(nPC_LP)','.pdf'];
      
      % Normilize Variance of SHOCK
      if norm_SHOCK == 1
            sdSHOCK          = nanstd(SHOCK);
            SHOCK            = SHOCK/sdSHOCK;
      end
      
      %Initializating the loop
      for kk = 1:size(dep_var,2)
            % Define inputs for local_projection
            varnamekk                   = varlist{kk};
            disp(['  - Projecting ',varnamekk])
            fprintf('\n')
            depvarkk                    = dep_var(:,kk);
            [~, loc_start(kk,is), loc_end]     = truncate_data([depvarkk SHOCK PC_LP]);
            % Starting the Sample after loc_start_exogenous
            lockk = find(loc_start_exogenous == Time);
            if loc_start(kk,is) < lockk
                  loc_start(kk,is) = lockk;
            end
            depvarkk                    = depvarkk(loc_start:loc_end);
            SHOCKkk                     = SHOCK(loc_start:loc_end);
            pckk                        = PC_LP(loc_start:loc_end,:);

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
      plot_IRF(varlist,IRF_low,IRF_low2,IRF_up,IRF_up2,IRF(:,:,is),H,printIRFs,IRFname); %change this function
      % Plot VD
      plot_IRF(varlist,VDkk,VDkk,VDkk,VDkk,VDkk,H,printVD,VDname)
      close
end

% Technical Parameters
tech_info_table;

return


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

%Compute average spectral density, D1, around the peak  and average
%spectral density around the trough, D2
lpeak_low    = 24; %should be adjusted with steps and IRF horizon
lpeak_up     = 26;
ltrough_low  = 58;
ltrough_up   = 60;
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
