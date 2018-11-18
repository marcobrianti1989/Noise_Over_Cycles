clear
close all

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
lag                 = 2;  %number of lags of control variables (other structural shocks)
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
pause(0.5)
close


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

% Standardize Ztilde to get one std dev shock
Ztilde  = Ztilde/std(Ztilde);
Ztilde  = [nan(loc_start-1,1); Ztilde; nan(size(dataset,1)-loc_end,1)];

% Matrix of dependen variables - All the variables are in log levels
control_pop = 0; % Divide GDP, Cons, Hours, Investment over population
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
    loc_start                   = loc_start + lags;
    depvarkk                    = depvarkk(loc_start:loc_end);
    Ztildekk                    = Ztilde(loc_start:loc_end);
    % Run local_projection
    [IR{kk},res{kk},Rsquared{kk},BL{kk},tuple{kk},VarY{kk}] = ...
        local_projection(depvarkk,pckk,Ztildekk,lags,H);
    [B, Yhat, res]  = quick_ols(Y,Ztildekk);
    XX             = [XX resEBP(1+nlags:end) resJLN(1+nlags:end)];%
    LM             = fitlm(XX,Y,'linear')
end

