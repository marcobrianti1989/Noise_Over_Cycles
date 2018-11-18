%*************************************************************************%
% Main
%
% NOTE describe variables (especially SHOCKS) in dataset
%
% last change 8/17/2018
%
% Code by Brianti, Marco e Cormun, Vito
%*************************************************************************%

clear
close all

%Read main dataset
filename                    = 'main_file';
sheet                       = 'Sheet1';
range                       = 'B1:BJ300';
do_truncation               = 0; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data2(filename, sheet, range, do_truncation);
dataset                     = real(dataset);
% numberTFP                   = strmatch('DTFP_UTIL', var_names);
%Assess names to each variable
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

%Read dataset_PC for PC analysis
filename_PC                                = 'Dataset_test_PC';
sheet_PC                                   = 'Quarterly';
range_PC                                   = 'B2:DA300';
do_truncation_PC                           = 1; %Do truncate data.
[dataset_PC, var_names_PC]                 = read_data2(filename_PC, sheet_PC, range_PC, do_truncation_PC);
dataset_PC                                 = real(dataset_PC);
date_start_PC                              = dataset_PC(1,1);
dataset_PC                                 = dataset_PC(:,2:end); %Removing time before PC analysis
Zscore                                     = 1; %Standardize data before taking PC
PC                                         = get_principal_components(dataset_PC,Zscore);
pc                                         = nan(size(dataset,1),size(dataset_PC,2));
loc_time_PC                                = find(Time == date_start_PC);
pc(loc_time_PC:loc_time_PC+size(PC,1)-1,:) = PC;

%*************************************************************************%
%                                                                         %
%          1st stage - Deriving noise shocks                              %
%                                                                         %
%*************************************************************************%

%Building Zt
%Step 1 - Getting the forecasted growth rates
%Real GDP
Delta_RGDP_t        = log(RGDP5_SPF) - log(RGDP1_SPF);
Delta_RDGP_t1       = log(RGDP6_SPF) - log(RGDP2_SPF);
%Industrial Production
Delta_INDPROD_t     = log(dataset(:,22)) - log(dataset(:,20));
Delta_INDPROD_t1    = log(dataset(:,23)) - log(dataset(:,21));
%Investment is the sum between residential and non residential investment
Delta_RINV_t        = log(dataset(:,14) + dataset(:,18)) - log(dataset(:,12) + dataset(:,16));
Delta_RINV_t1       = log(dataset(:,15) + dataset(:,19)) - log(dataset(:,13) + dataset(:,17));
%Step 2 - Revision in forecast growth rates
Z1                  = [NaN; Delta_RGDP_t(2:end) - Delta_RDGP_t1(1:end-1)];
Z2                  = [NaN; Delta_INDPROD_t(2:end) - Delta_INDPROD_t1(1:end-1)];
Z3                  = [NaN; Delta_RINV_t(2:end) - Delta_RINV_t1(1:end-1)];
ZZ                  = Z1; %Select GDP growth

%Technical values to build Ztilde
lag_tfp             = 4; %number of lags of TFP - cannot be zero since 1 include current TFP
lead_tfp            = 16; %number of leads of TFP
lag                 = 0;  %number of lags of control variables (other structural shocks)
mpc                 = 2; %max number of principal components
threshold           = -1/eps; %Remove all the NaN values

%Runniong OLS to obtain Ztilde
T                 = size(ZZ,1);
const             = ones(T,1);
X                 = const;

%Structural shocks from Ramey narrative approach, Hamilton, Romer and
%Romer, Military government spending...
controls                    = [MUNI1Y,PDVMILY,HAMILTON3YP,RESID08,TAXNARRATIVE];
trend                       = 1:1:length(X); %Control for the time trend
X                           = [X, trend', controls];
[data, loc_start, loc_end]  = truncate_data([ZZ X]);
loc_start                   = loc_start + lag;
ZZ                          = data(lag+1:end,1);
X                           = data(lag+1:end,2:end);
DTFP                        = [NaN; diff(TFP)];

%Control for TFP
for i = 1:lag_tfp %Add lags of TFP - When i = 1 TFP is contemporaneous
      X(:,end+1)  = DTFP(loc_start+1-i:loc_end-i+1);
end
for i = 1:lead_tfp %Add leads of TFP
      X(:,end+1)  = DTFP(loc_start+i:loc_end+i);
end
for l = 1:lag %Add lags of controls
      X           = [X controls(loc_start-l:loc_end-l,:) ...
            pc(loc_start-l:loc_end-l,1:mpc)];
end
Y                 = ZZ;
[B,zhat,Ztilde]   = quick_ols(Y,X);

%Show the graph of Ztilde - Figure(1)
plot1 = 1; % if plot = 1, figure will be displayed
plot_Ztilde(Ztilde,Time,NBERDates,loc_start,loc_end,plot1)

%Correlation with Barsky and Sims 2011

news = BarskySims_News(loc_start:loc_end);
corr_ztilde_news = corr(Ztilde,news);


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

% Create Var List
% varlist          = {'TFP','RealGDP', 'RealCons',...
%       'UnempRate','RealWage','Hours','CPIInflation',...
%       'RealInvestment','SP500','OilPrice','GZSpread','FFR',...
% 'Vix','VXO','Inventories','LaborProductivity','Spread'};
varlist          = {'RealGDP', 'RealCons','UnempRate','Hours','RealInvestment',...
      'RealInventories','RealProfitsaT','RealSales',...
      'CPIInflation','PriceCPE','CPIDurables',... %All the nominal variables should be last
      'CPINonDurables'};
numberCPI        = strmatch('CPIInflation', varlist);
numberCPE        = strmatch('PriceCPE', varlist);
numberCPID       = strmatch('CPIDurables', varlist);
numberCPIND      = strmatch('CPINonDurables', varlist);
numberCPIS       = strmatch('CPIServices', varlist);
numberGDP        = strmatch('RealGDP', varlist);
numberC          = strmatch('RealCons', varlist);
numberHours      = strmatch('Hours', varlist);
numberInv        = strmatch('RealInvestment', varlist);
numberProf       = strmatch('RealProfitsaT', varlist);
numberInvent     = strmatch('RealInventories', varlist);

%numberInflation  = strmatch('Inflation', varlist);
lags             = 3;
H                = 20; %irfs horizon
mpc              = 3; %max number of principal components

% Standardize Ztilde to get one std dev shock
Ztilde  = Ztilde/std(Ztilde);
Ztilde  = [nan(loc_start-1,1); Ztilde; nan(size(dataset,1)-loc_end,1)];

% Matrix of dependen variables - All the variables are in log levels
control_pop = 1; % Divide GDP, Cons, Hours, Investment over population
for i = 1:length(varlist)
      dep_var(:,i) = eval(varlist{i});
      if control_pop == 1
            if i == numberGDP || i == numberC || i == numberHours || i == numberInv || i == numberInvent || i == numberInvent
                  dep_var(:,i) = dep_var(:,i) - Population;
            end
      end
end

% Set up year on year inflation
use_Inflation = 1;
if use_Inflation == 1
      for ii = 1:length(dep_var)-4
            CPI(ii)       = dep_var(ii+4,numberCPI) - dep_var(ii,numberCPI);
            CPE(ii)       = dep_var(ii+4,numberCPE) - dep_var(ii,numberCPE);
            CPID(ii)      = dep_var(ii+4,numberCPID) - dep_var(ii,numberCPID);
            CPIND(ii)     = dep_var(ii+4,numberCPIND) - dep_var(ii,numberCPIND);
            %CPIS(ii)      = dep_var(ii+4,numberCPIS) - dep_var(ii,numberCPIS);
      end
      CPI     = [NaN; NaN; NaN; NaN; CPI'];
      CPE     = [NaN; NaN; NaN; NaN; CPE'];
      CPID    = [NaN; NaN; NaN; NaN; CPID'];
      CPIND   = [NaN; NaN; NaN; NaN; CPIND'];
      %CPIS    = [NaN; NaN; NaN; NaN; CPIS'];
      loc     = min([numberCPI numberCPE numberCPID numberCPIND numberCPIS]);
      dep_var = [dep_var(:,1:loc-1) CPI CPE CPID CPIND]; %CPIS];
end

% Set up the typology of transformation
logdifferences = 0;
if logdifferences == 1
      dep_var = [nan(1,size(dep_var,2)); diff(dep_var)];
end

for kk = 1:size(dep_var,2)
      % Define inputs for local_projection
      depvarkk                    = dep_var(:,kk);
      [~, loc_start, loc_end]     = truncate_data([depvarkk Ztilde pc]);
      loc_start                   = loc_start; %+ lags;
      depvarkk                    = depvarkk(loc_start:loc_end);
      Ztildekk                    = Ztilde(loc_start:loc_end);
      pckk                        = pc(loc_start:loc_end,1:mpc);
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
      nsimul         = 200;
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
n_row    = 3; % how many row in the figure
unique   = 1; % if unique = 1 plot IRFs together, if = 1 plot each IRF separately
plot_IRF_lp_unconditional(varlist,IRF_low,IRF_low2,IRF_up,IRF_up2,IRF,H,plot2,n_row,unique)

% Print figure authomatically if "export_figure1 = 1"
if plot2 == 1
      export_fig2 = 0; % if export_fig1 = 1, figure will be saved
      export_fig_IRF_lp_unconditional(export_fig2)
end

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
















