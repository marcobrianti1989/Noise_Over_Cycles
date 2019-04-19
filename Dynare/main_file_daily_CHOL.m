%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Marco Brianti, PhD Candidate, Boston College, Department of Economics, Dec 4, 2018
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            READING DATA                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reading Data
filename                    = 'main_file_daily';
sheet                       = 'daily';
range                       = 'B1:V1535';
do_truncation               = 0; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data(filename, sheet, range, do_truncation);
tf                          = isreal(dataset);
if tf == 0
      fprintf('\n')
      warning('Dataset has complex variables in it.')
      fprintf('\n')
end
dataset                     = real(dataset);
% time_start                  = dataset(1,1);
% time_end                    = dataset(end,1);
% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

% Read Days
filename                    = 'main_file_daily';
sheet                       = 'daily';
range                       = 'A7:A1535';
[~, DAYS]                   = xlsread(filename,sheet,range);


%showCDS(CDSITA2014,DAYS,Time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          CREATING SHOCKS                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loading Keydates
keydates_script;

% Create Instrument
start_later = 1;
dummyITA = create_instrument(Time,keydates,DAYS,start_later,CDSITA,CDSITA2014);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          STRUCTURAL VAR                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the system
system_names    = {'dummyITA','CDSITA2014','CDSBANKS2014','IVFTSE','FTSE'};
loc_dummy       = find(strcmp(system_names,'dummyITA') == 1);
loc_CDSITA      = find(strcmp(system_names,'CDSITA') == 1);

% Build System and eventually filter it
for i = 1:length(system_names)
      system(:,i) = eval(system_names{i});
end

% Truncate data and create the system variables
%start_later    = 1;
start_date     = '24SEP2014';
loc_start_date = find(strcmp(DAYS,start_date) == 1);
if start_later == 1
      system = system(loc_start_date:end,:);
      Time   = Time(loc_start_date:end,:);
      %DAYS   = DAYS{loc_start_date:end};
      disp('Data are starting in September 2014.')
      fprintf('\n')
else
      disp('Data are starting in January 2013.')
      fprintf('\n')
end

% Detrend system
which_trend = 'demean'; %BP, HP, lin, diff, none, demean
meansystem  = mean(system,1);
system      = detrend_func(system,which_trend);

% Tests for lags
max_lags     = 40;
[AIC,BIC,HQ] = aic_bic_hq(system,max_lags);

% Cholesky decomposition
nlags           = 2;
disp(['Number of lags is ',num2str(nlags)])
fprintf('\n')
[A,B,res,sigma] = sr_var(system, nlags);

% Create dataset from bootstrap
nburn             = 0;
nsimul            = 20;
which_correction  = 'none';
blocksize         = 4;
[beta_tilde, data_boot2, beta_tilde_star,nonstationarities] ...
      = bootstrap_with_kilian(B,nburn,res,nsimul,which_correction,blocksize);

% Get "bootstrapped A" nsimul times
for i_simul=1:nsimul
      % Cholesky decomposition
      [A_boot(:,:,i_simul),B_boot(:,:,i_simul),~,~] = ...
            sr_var(data_boot2(:,:,i_simul), nlags);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      IMPULSE RESPONSE FUNCTIONS                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate IRFs with upper and lower bounds
sig1                       = 0.05;
sig2                       = 0.025;
H                          = 22;
normIRFs                   = 1;
[IRFs, ub1, lb1, ub2, lb2] = genIRFs(A,A_boot,B,B_boot,H,sig1,sig2,normIRFs);

IRFs = IRFs(2:end,:,:);
ub1 = ub1(2:end,:,:);
lb1 = lb1(2:end,:,:);
ub2 = ub2(2:end,:,:);
lb2 = lb2(2:end,:,:);

% Create and Printing figures for IRFs
base_path         = pwd;
which_ID          = 'IRFs_2LAG';
print_figs        = 'no';
use_current_time  = 1; % (don't) save the time
which_shocks      = [loc_dummy]; %[Uposition];
shocknames        = {'Uncertainty Shock'};
unique            = 1;
system_names      = {'CDSITA2014','CDSBANKS2014','IVFTSE','FTSE'};
load IRFIVLP_2LAGS_SEP2014
load IRFIVSVAR_2LAGS_SEP2014
if unique == 1
      plot_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,H,which_shocks,shocknames,...
            system_names,which_ID,print_figs,use_current_time,base_path,0,0)
else
      plot_IRFs_2CIs_multifigures(IRFs,ub1,lb1,ub2,lb2,H,which_shocks,shocknames,...
            system_names,which_ID,print_figs,use_current_time,base_path)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       VARIANCE DECOMPOSITION                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create and Printing figures for Variance decomposition
which_ID          = 'vardec_2LAG';
print_figs        = 'no';
sig1VD            = 0.05;
sig2VD            = 0.025;
normVD            = 0;
[vardec, ub1_vardec, lb1_vardec, ub2_vardec, lb2_vardec] = ...
      gen_vardec_boot(A,A_boot,A,A_boot,B,B_boot,H,sig1VD,sig2VD,normVD);
% Plotting VD
plot_IRFs_2CIs(vardec(2:end,:,:),ub1_vardec(2:end,:,:),lb1_vardec(2:end,:,:),...
      ub2_vardec(2:end,:,:),lb2_vardec(2:end,:,:),H,which_shocks,shocknames,...
      system_names,which_ID,print_figs,use_current_time,base_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       STRUCTURAL SHOCKS                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get Structural Shocks
ss  = (inv(A)*res')';
ssU = ss(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     HISTORICAL DECOMPOSITION                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Historical Decomposition
normHD          = 0;
[HD,diffHD]     = historical_decomposition(A,B,ss,normHD);
HD              = HD(2:end,:,:);
system_namesHD  = {'CDSITA2014','CDSBANKS2014','IVFTSE','FTSE'};
systemHD        = system(:,2:end);
meansystemHD    = meansystem(2:end);

% Plotting Figures
which_ID          = 'HD_2LAG';
print_figs        = 'yes';
unique            = 0;
year_to_start     = 0;
which_shocksHD    = loc_dummy;
plot_historical_decomposition(Time,HD,systemHD,meansystemHD,DAYS,nlags,...
      which_shocks,shocknames,system_namesHD,print_figs,...
      use_current_time,base_path,unique,which_ID,year_to_start)


asd


