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
     
% Define inputs for local_projection
[~, loc_start, loc_end]     = truncate_data([RGDPG ZZ Delta_RGDP_t]);
loc_start                   = loc_start + 0;
loc_end                     = loc_end - 0;
RGDPG                       = RGDPG(loc_start:loc_end);
ZZ                          = ZZ(loc_start:loc_end);
Delta_RGDP_t                = Delta_RGDP_t(loc_start:loc_end);
FE                          = [NaN; NaN; NaN; NaN; RGDPG(1+4:end) - Delta_RGDP_t(1:end-4)];

% Run OLS of forecast error on the forecast revision
[B, Yhat, res]  = quick_ols(FE(1+4:end),ZZ(1:end-4));
LM              = fitlm(ZZ(1:end-4),FE(1+4:end),'linear')
% scatter(RGDPG(1+4:end),Delta_RGDP_t(1:end-4))
scatter(FE(1+4:end),ZZ(1:end-4))

