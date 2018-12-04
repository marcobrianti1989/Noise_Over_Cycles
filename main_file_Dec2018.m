%*************************************************************************%
% Main
%
% NOTE describe variables (especially SHOCKS) in dataset
%
% last change 11/30/2018
%
% Code by Brianti, Marco e Cormun, Vito
%*************************************************************************%

clear
close all

% Read main dataset
filename                    = 'main_file';
sheet                       = 'Sheet1';
range                       = 'B1:BU300';
do_truncation               = 0; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data2(filename, sheet, range, do_truncation);
nNaN                        = 20; % adding some NaN at the end to have space for leads
dataset                     = [dataset; NaN(nNaN,size(dataset,2))];
% Assess name to each variable
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

%*************************************************************************%
%                                                                         %
%          1st stage - Deriving sentimentcks                              %
%                                                                         %
%*************************************************************************%

%Building Zt
%Step 1 - Getting the forecasted growth rates
Delta_RGDP_t        = log(RGDP5_SPF) - log(RGDP1_SPF);
Delta_RDGP_t1       = log(RGDP6_SPF) - log(RGDP2_SPF);
%Step 2 - Revision in forecast growth rates
Z                   = [NaN; Delta_RGDP_t(2:end) - Delta_RDGP_t1(1:end-1)];

% Building the Forecast Error of
%RGDP_known_t        = log(RGDP1_SPF(2:end)); % from 2 as it is RGDP at t, infoset (t+1)
%RGDP_forec          = log(RGDP5_SPF(1:end-1)); % infoset at t, Forecast t+3
FE                  = [NaN(3,1); log(RGDP1_SPF(1+4:end)) - log(RGDP5_SPF(1:end-4)); NaN];

%Technical values to build Ztilde
lags                 = 4; %number of lags of TFP - cannot be zero since 1 include current TFP
leads                = 16; %number of leads of TFP

%Structural shocks from Ramey narrative approach, Hamilton, Romer and Romer, Military government spending...
[TFP_trunc, trunc1, trunc2] = truncate_data(TFP);
TFPBP                       = bpass(TFP_trunc,4,32);
TFPBP                       = [TFPBP; NaN(length(TFP) - length(TFPBP),1)];
PC                          = [PC1 PC2 PC3];
X_contemporaneous           = [TFPBP MUNI1Y PDVMILY HAMILTON3YP RESID08 TAXNARRATIVE];
X_lag                       = [TFPBP PC MUNI1Y PDVMILY HAMILTON3YP RESID08 TAXNARRATIVE];
X_lead                      = TFPBP;
Y                           = Z;

% Control Regression
[~, ~, Ztilde,regressor] = lead_lag_matrix_regression(Y,X_lead,leads,X_lag,lags,...
      X_contemporaneous);

%Show the graph of Ztilde - Figure(1)
plot1 = 1; % if plot = 1, figure will be displayed
plot_Ztilde(Ztilde,Time(1+lags:end-leads),NBERDates(1+lags:end-leads),plot1)

% Print figure authomatically if "export_figure1 = 1"
if plot1 == 1
      export_fig1 = 0; % if export_fig1 = 1, figure will be saved
      export_fig_Ztilde(export_fig1)
end

%*************************************************************************%
%                                                                         %
%          2nd stage - Smooth Transition Local Projections                %
%                                                                         %
%*************************************************************************%

% Preparing Data Dependent Variables

% Create Inflation from Price Indexes
ninfl = 4;
PCEInflation    = create_inflation(PriceCPE,ninfl);
CPIInflation    = create_inflation(CPIInflation,ninfl);
CPIDurables     = create_inflation(CPIDurables,ninfl);
CPINonDurables  = create_inflation(CPINonDurables,ninfl);
CPIServices     = create_inflation(CPIServices,ninfl);

% Per capita adjustment
control_pop = 0; % Divide GDP, Cons, Hours, Investment over population
if control_pop == 1
      RealGDP                 = RealGDP - Population;
      RealCons                = RealCons - Population;
      RealInvestment          = RealInvestment - Population;
      Hours                   = Hours + Employment - Population;
      RealInventories         = RealInventories - Population;
      RealSales               = RealSales - Population;
end

% Create Var List
SP500            = SP500 - GDPDefl;
varlist          = {'RealGDP', 'RealCons','SP500',...
      'Hours','RealInvestment','RealInventories',...
      'TFP','UnempRate','RealSales',...
      'FE','CPIInflation','PCEInflation'};

% Matrix of dependen variables - All the variables are in log levels
for i = 1:length(varlist)
      if strcmp(varlist{i},'Ztilde') == 1
            dep_var(:,i) = [NaN(lags,1); eval(varlist{i}); NaN(leads,1)];
      else
            dep_var(:,i) = eval(varlist{i});
      end
end

% Set up the typology of transformation
logdifferences = 0;
if logdifferences == 1
      dep_var = [nan(1,size(dep_var,2)); diff(dep_var)];
end

% Align the timing - It is important!
dep_var          = dep_var(1+lags:end-leads,:);
PC               = PC(1+lags:end-leads,:);

% Technical Parameters
H                = 20; %irfs horizon
lags             = 4;
HPfilter         = 0;
BPfilter         = 0;
sdZtilde         = nanstd(Ztilde);
Ztilde           = Ztilde/sdZtilde;

for kk = 1:size(dep_var,2)
      % Define inputs for local_projection
      depvarkk                    = dep_var(:,kk);
      [~, loc_start, loc_end]     = truncate_data([depvarkk Ztilde PC]);
      depvarkk                    = depvarkk(loc_start:loc_end);
      if HPfilter == 1
            if strcmp('FE',varlist{kk}) == 1
                  disp('FE is not HP filtered')
            else
                  [~, depvarkk]         = hpfilter(depvarkk,1600);
            end
      end
      if BPfilter == 1
            if strcmp('FE',varlist{kk}) == 1
                  disp('FE is not BP filtered')
            else
                  depvarkk              = bpass(depvarkk,4,32);
            end
      end
      Ztildekk                    = Ztilde(loc_start:loc_end);
      pckk                        = PC(loc_start:loc_end,:);
      % Run local_projection
      [IR{kk},res{kk},Rsquared{kk},BL{kk},tuple{kk},VarY{kk}] = ...
            local_projection(depvarkk,pckk,Ztildekk,lags,H);
      if logdifferences == 0
            IRF(kk,:) = IR{kk};
      else
            IRF(kk,:) = cumsum(IR{kk});
      end
      % Build a table for the Variance Explained by Ztilde - Following  Stock,
      % Watson (2018) - The Economic Journal, page 928 Eq. (15)
      VarY_ih = VarY{kk};
      for ih = 1:H
            VarYY    = VarY_ih(ih);
            VarExplained(kk,ih) = sum(IRF(kk,1:ih).^2)/VarYY;
      end
      % Initiate bootstrap
      nsimul         = 500;
      tuplekk        = tuple{kk};
      for hh = 1:H
            tuplekkhh = tuplekk{hh}; % Fix a specific horizon
            Y                             = tuplekkhh(:,1);
            X                             = tuplekkhh(:,2:end);
            XControl                      = tuplekkhh(:,3:end);
            [Yboot, Xboot]                = bb_bootstrap_LP(Y,X,nsimul,lags);
            [YbootC, XbootC]              = bb_bootstrap_LP(Y,XControl,nsimul,lags);
            for isimul = 1:nsimul
                  B                       = Xboot(:,:,isimul)'*Xboot(:,:,isimul)\...
                        (Xboot(:,:,isimul)'*Yboot(:,isimul));
                  BC                      = XbootC(:,:,isimul)'*XbootC(:,:,isimul)\...
                        (XbootC(:,:,isimul)'*YbootC(:,isimul));
                  IRF_boot(kk,hh,isimul)  = B(1);
                  VarYBoot(kk,hh,isimul)  = var(YbootC(:,isimul) - XbootC(:,:,isimul)*BC);
            end
      end
end

% Select upper and lower bands
for kk = 1:size(dep_var,2)
      IRF_bootkk = IRF_boot(kk,:,:);
      VarYbootkk = VarYBoot(kk,:,:);
      if logdifferences == 0
            IRF_boot(kk,:,:)  = IRF_bootkk;
            VarY_boot(kk,:,:) = VarYbootkk;
      else
            IRF_boot(kk,:,:)  = cumsum(IRF_bootkk,2);
            VarY_boot(kk,:,:) = cumsum(VarYbootkk,2);
      end
end
IRF_boot         = sort(IRF_boot,3);
VarY_boot        = sort(VarY_boot,3);
sig              = 0.05;
sig2             = 0.16;
up_bound         = floor(nsimul*sig); % the upper percentile of bootstrapped responses for CI
up_bound2        = floor(nsimul*sig2); % the upper percentile of bootstrapped responses for CI
low_bound        = ceil(nsimul*(1-sig)); % the lower percentile of bootstrapped responses for CI
low_bound2       = ceil(nsimul*(1-sig2)); % the lower percentile of bootstrapped responses for CI
IRF_up           = IRF_boot(:,:,up_bound);
VarY_up          = VarY_boot(:,:,up_bound);
IRF_up2          = IRF_boot(:,:,up_bound2);
VarY_up2         = VarY_boot(:,:,up_bound2);
IRF_low          = IRF_boot(:,:,low_bound);
VarY_low         = VarY_boot(:,:,low_bound);
IRF_low2         = IRF_boot(:,:,low_bound2);
VarY_low2        = VarY_boot(:,:,low_bound2);

% Confidence Intervals for Variance Explained
for kk = 1:size(dep_var,2)
      VarYup   = VarY_up(kk,:);
      VarYup2  = VarY_up2(kk,:);
      VarYlow  = VarY_low(kk,:);
      VarYlow2 = VarY_low2(kk,:);
      for ih = 1:H
            VarYYup   = VarYup(ih);
            VarYYup2  = VarYup2(ih);
            VarYYlow  = VarYlow(ih);
            VarYYlow2 = VarYlow2(ih);
            VarExplainedup(kk,ih)   = sum(IRF_up(kk,1:ih).^2)/VarYYup;
            VarExplainedlow(kk,ih)  = sum(IRF_low(kk,1:ih).^2)/VarYYlow;
            VarExplainedup2(kk,ih)  = sum(IRF_up2(kk,1:ih).^2)/VarYYup2;
            VarExplainedlow2(kk,ih) = sum(IRF_low2(kk,1:ih).^2)/VarYYlow2;
      end
end

%Show the graph of IRF - Figure(2)
plot2    = 1; % if plot2 = 1, figure will be displayed
n_row    = 2; % how many row in the figure
unique   = 1; % if unique = 1 plot IRFs together, if = 1 plot each IRF separately
plot_IRF_lp_unconditional(varlist,100.*IRF_low,100.*IRF_low2,100.*IRF_up,100.*IRF_up2,100.*IRF,H,plot2,n_row,unique)

%Print figure authomatically if "export_figure1 = 1"
if plot2 == 1
      export_fig2 = 0; % if export_fig1 = 1, figure will be saved
      export_fig_IRF_lp_unconditional(export_fig2)
end

asdf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variance Decomposition - Need to be checked!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Show the variance Explained - Figure(3)
plot3    = 1; % if plot2 = 1, figure will be displayed
n_row    = 3; % how many row in the figure
unique   = 1; % if unique = 1 plot IRFs together, if = 1 plot each IRF separately
plot_IRF_lp_unconditional(varlist,VarExplained,VarExplained,VarExplained,...
      VarExplained,VarExplained,H,plot3,n_row,unique)

% Print figure authomatically if "export_figure1 = 1"
if plot3 == 1
      export_fig3 = 0; % if export_fig1 = 1, figure will be saved
      export_fig_IRF_lp_unconditional(export_fig3)
end
















