%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Marco Brianti and Vito Cormun,
%  PhD Candidates, Boston College,
%  Department of Economics, February 27, 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

% Technical parameters
nlags                      = 8;      %lags in the VAR system
control_pop                = 0;       % Divide key macro variables over population
which_trend                = 'quad';  %'BP', 'HP', 'lin', 'quad', 'diff', 'none', 'demean': detrending the variables before adding them in the VAR.
which_boot                 = 'none';  % Either 'none' or 'blocks'
blocksize                  = 4;       % if which_boot = 'blocks', then decide the block of residuals
nsimul                     = 20;     % Number of bootstrap simulations
nburn                      = 0;       % Number of observations to burn during each bootstrap
sig1                       = 0.16;    % Tighter Confidence Interval
sig2                       = 0.05;    % Looser Confidence Interval
H                          = 25;      % Horizon of IRFs
print_figs                 = 'yes';    % If you want to pring fig 'yes', otherwise 'no'
use_current_time           = 1;       % In order to avoid overwriting
shocknames                 = {'Sentiment Shock'};                     % Name of the Shock(s)
system_names               = {'ZTILDE','RealGDP','RealInvestment','HoursAll',...
      'RealInventories','RealCons','GDPDefl'};             % Variables in the VAR  %'RealGDP','RealCons','RealInvestment','HoursPerPerson'};
pos_ZTILDE                 = find(strcmp('ZTILDE',system_names)==1);  % Position of Ztilde
pos_uTFP                   = find(strcmp('UnantTFPshock',system_names)==1); % Position of unexpected productivity shock
which_shocks               = pos_ZTILDE;                              % Cholesky, Position of the Shock
system_names_graph         = {'Real GDP','Real Investment','Total Hours',...
      'Real Inventories','Real Consumption','GDP Deflator'};  
% Read main dataset
filename                    = 'main_file';
sheet                       = 'Sheet1';
range                       = 'B1:DX300';
do_truncation               = 0; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data2(filename, sheet, range, do_truncation);
tf                          = isreal(dataset);
if tf == 0
      warning('Dataset has complex variables in it.')
end
dataset                     = real(dataset);
time_start                  = dataset(1,1);
time_end                    = dataset(end,1);
% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

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

subsystem_names{1} = system_names{1};
for iss = 1:length(system_names)-1
      subsystem_names{2} = system_names{iss+1};
      % Create VAR System and Truncate NaN
      for i = 1:length(subsystem_names)
            subsystem(:,i) = eval(subsystem_names{i});
      end
      [subsystem, truncation_point, truncation_point2] = truncate_data(subsystem);
      
      % Detrend Variables
      subsystem       = detrend_func(subsystem,which_trend);
      
      % Tests for lags
      max_lags     = 8;
      [AIC(iss),BIC(iss),HQ(iss)] = aic_bic_hq(subsystem,max_lags);
      
      % Cholesky decomposition
      [A(:,:,iss),B(:,:,iss),res,sigma(:,:,iss)] = sr_var(subsystem, nlags);
      RES{iss}                                   = res;
      
      % Get Structural Shocks
      ss          = (inv(A(:,:,iss))*res')';
      SS{iss}     = ss;
      
      % Create dataset from bootstrap
      [~, data_boot2,~,~] ...
            = bootstrap_with_kilian(B(:,:,iss),nburn,res,nsimul,which_boot,blocksize);
      
      % Get A and B with Cholesky decomposition for each simulated series
      for i_simul=1:nsimul
            [A_boot(:,:,i_simul,iss),B_boot(:,:,i_simul,iss),~,~] = ...
                  sr_var(data_boot2(:,:,i_simul),nlags);
      end
      
      % Generate IRFs with upper and lower bounds
      [IRFs(:,:,:,iss), ub1(:,:,:,iss), lb1(:,:,:,iss), ub2(:,:,:,iss), lb2(:,:,:,iss)] = ...
            genIRFs(A(:,:,iss),A_boot(:,:,:,iss),B(:,:,iss),B_boot(:,:,:,iss),H,sig1,sig2);
      
      clear subsystem data_boot2 res ss
end

% Reshape IRFs to use "plot_IRFs_2CIs" function below
IRFS         = IRFs(2,:,1,:);
IRFS         = 100*squeeze(IRFS)';
UB1          = ub1(2,:,1,:);
UB1          = 100*squeeze(UB1)';
UB2          = ub2(2,:,1,:);
UB2          = 100*squeeze(UB2)';
LB1          = lb1(2,:,1,:);
LB1          = 100*squeeze(LB1)';
LB2          = lb2(2,:,1,:);
LB2          = 100*squeeze(LB2)';

% Create and Printing figures
for iend = 1:length(system_names)-1
      endogenous_var_names{iend} = system_names{iend+1};
end
base_path         = pwd;
which_ID          = 'CHOL';
which_shocks      = 1;
plot_IRFs_2CIs(IRFS,UB1,LB1,UB2,LB2,H,which_shocks,shocknames,...
      system_names_graph,which_ID,print_figs,use_current_time,base_path)

tech_info_table_Chol;










