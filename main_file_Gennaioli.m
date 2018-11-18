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

% Building ZZ
% Step 1 - Getting the forecasted GDP growth rates
Delta_RGDP_t                = log(RGDP5_SPF) - log(RGDP1_SPF);
Delta_RDGP_t1               = log(RGDP6_SPF) - log(RGDP2_SPF);
% Step 2 - Revision in forecast growth rates
ZZ                          = [NaN; Delta_RGDP_t(2:end) - Delta_RDGP_t1(1:end-1)];

% Compute realized real gdp growth
RGDPG                       = [NaN; NaN; NaN; NaN; RealGDP(1+4:end) - RealGDP(1:end-4)];

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
[data, loc_start, loc_end]  = truncate_data([ZZ X RGDPG]);
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
      X           = [X controls(loc_start-l:loc_end-l,:) pc(loc_start-l:loc_end-l,1:mpc)];
end
Y                 = ZZ;
[B,zhat,Ztilde]   = quick_ols(Y,X);

% Align RGDPG
RGDPG             = RGDPG(loc_start:loc_end);
Delta_RGDP_t      = Delta_RGDP_t(loc_start:loc_end);
FE                = [NaN; NaN; NaN; NaN; RGDPG(1+4:end) - Delta_RGDP_t(2:end-3)];

% Run OLS of forecast error on the forecast revision
[B, Yhat, res]  = quick_ols(FE(1+4:end),[ones(length(ZZ)-4,1) Ztilde(2:end-3)]);
LM              = fitlm(Ztilde(2:end-3),abs(FE(1+4:end)),'linear')
% scatter(RGDPG(1+4:end),Delta_RGDP_t(1:end-4))
% scatter(FE(1+4:end),ZZ(1:end-4))

figure
plot(6*Ztilde(2:end-3),'linewidth',2,'color','r')
hold on
plot(FE(1+4:end),'linewidth',2,'color','b')
plot(RGDPG(1+4:end),'linewidth',2,'color','k')
legend('Ztilde','FE','RGDPG')





