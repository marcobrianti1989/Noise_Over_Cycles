%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Marco Brianti, PhD Candidate, Boston College, Department of Economics, August 8, 2018
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

% Read main dataset
filename                    = 'main_file';
sheet                       = 'Sheet1';
range                       = 'B1:DP300';
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
control_pop = 0; % Divide GDP, Cons, Hours, Investment over population
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

%Building Zt - Forecast Revisions from SPF and Michigan Index
create_Z;
Z = Z1;

% Define the system1
system_names  = {'TFP','Z','RealGDP','RealCons','RealInvestment',...
      'HoursAll'};

for i = 1:length(system_names)
      system(:,i) = eval(system_names{i});
end
Zpos    = find(strcmp('Z', system_names));
[system, truncation_point, truncation_point2] = truncate_data(system);

% Detrend Variables
% which_trend = 'quad';
% system = detrend_func(system,which_trend);

% Tests for lags
max_lags     = 4;
[AIC,BIC,HQ] = aic_bic_hq(system,max_lags);

% Cholesky decomposition
nlags           = 1;
[A,B,res,sigma] = sr_var(system, nlags);

% Get Structural Shocks
ss = (inv(A)*res')';

% Create dataset from bootstrap
nburn             = 0;
nsimul            = 500;
which_correction  = 'none';
blocksize         = 4;
[beta_tilde, data_boot2, beta_tilde_star,nonstationarities] ...
      = bootstrap_with_kilian(B,nburn,res,nsimul,which_correction,blocksize);

% Get "bootstrapped A" nsimul times
for i_simul=1:nsimul
      % Cholesky decomposition
      [A_boot(:,:,i_simul),B_boot(:,:,i_simul),~,~] = sr_var(data_boot2(:,:,i_simul), nlags);
end

% Generate IRFs with upper and lower bounds
sig1                       = 0.05;
sig2                       = 0.025;
H                          = 40;
[IRFs, ub1, lb1, ub2, lb2] = genIRFs(A,A_boot,B,B_boot,H,sig1,sig2);

% Create and Printing figures
base_path         = pwd;
which_ID          = 'chol_';
print_figs        = 'no';
use_current_time  = 1; % don't save the time
which_shocks      = [1]; %[Uposition];
shocknames        = {'Sentiment Shock'};


plot_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,H,which_shocks,shocknames,...
      system_names,which_ID,print_figs,use_current_time,base_path)
asd
% Get variance Decomposition
[IRF_vardec, ~, ~, ~, ~] = genIRFs(A,0,B,0,H,sig1,sig2);
m = linspace(1,H,H);
for im = 1:length(m)
      vardec(:,:,im) = gen_vardecomp(IRF_vardec,m(im),H);
end
vardec = vardec(TFPposition,4,:);
horz = linspace(0,H,H);
figure
hold on
plot(horz,vardec(TFPposition,:),'linewidth',2)
grid on
legend boxoff
xlabel('Horizon')
ylabel('Variance Explained')
title('Variance Explained Of Real GDP')

%save workspace_nicespecification_cons_inv_adjusted_Ulast







