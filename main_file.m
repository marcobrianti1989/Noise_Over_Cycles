%*************************************************************************%
% Main
%
% NOTE describe variables (especially SHOCKS) in dataset
%
% last change 6/21/2018
%
% Code by Brianti, Marco e Cormun, Vito
%*************************************************************************%

clear
close all

%Read main dataset
filename                    = 'main_file';
sheet                       = 'Sheet1';
range                       = 'B1:AU300';
do_truncation               = 0; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data2(filename, sheet, range, do_truncation);
dataset                     = real(dataset);
filter                      = 0; %Filter the data - Just HPfilter for now
numberTFP                   = strmatch('DTFP_UTIL', var_names);
%If filter = 0 check the excel file to set log-vediation, transformation = 5.
if filter == 1
    [dataset_filter, loc, loc2]               = truncate_data(dataset(:,numberTFP:end));
    loc2                                      = loc + loc2 - 1;
    switch 'hpfilter'
        case 'hpfilter'
            [~ , dataset_HP]                      = hpfilter(dataset_filter,1600);
            dataset_HP2(1:loc-1,:)                = NaN(loc-1,size(dataset_HP,2));
            dataset_HP2(loc:loc2,:)               = dataset_HP;
            dataset_HP2(loc2+1:size(dataset,1),:) = NaN(size(dataset,1)-loc2,size(dataset_HP,2));
            dataset                               = [dataset(:,1:numberTFP-1) dataset_HP2];
        case 'OneSidedHPfilter'
            dataset_HP                            = zeros(size(dataset_filter,1),size(dataset_filter,2));
            for iif = 1:length(dataset_filter)
                [~ , dataset_HPgen]             = hpfilter(dataset_filter(1:iif,:),1600);
                dataset_HP(iif,:)               = dataset_HPgen(end,:);
            end
            cutoff                                = 8;
            dataset_HP2(1:loc-1+cutoff,:)         = NaN(loc-1+cutoff,size(dataset_HP,2));
            dataset_HP2(loc+cutoff:loc2,:)        = dataset_HP(1+cutoff:end,:);
            dataset_HP2(loc2+1:size(dataset,1),:) = NaN(size(dataset,1)-loc2,size(dataset_HP,2));
            dataset                               = [dataset(:,1:numberTFP-1) dataset_HP2];
        case 'BandPassfilter' %Do not use it, still not working!
            dataset_bp                            = bpfilt(dataset_filter,6,32,1,0);
            dataset_bp2(1:loc-1,:)                = NaN(loc-1,size(dataset_bp,2));
            dataset_bp2(loc:loc2,:)               = dataset_bp(1:end,:);
            dataset_bp2(loc2+1:size(dataset,1),:) = NaN(size(dataset,1)-loc2,size(dataset_bp,2));
            dataset                               = [dataset(:,1:numberTFP-1) dataset_bp2];
        otherwise
            disp('Wrong filter selection')
    end
end
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
Z1                  = Delta_RGDP_t(2:end) - Delta_RDGP_t1(1:end-1);
Z2                  = Delta_INDPROD_t(2:end) - Delta_INDPROD_t1(1:end-1);
Z3                  = Delta_RINV_t(2:end) - Delta_RINV_t1(1:end-1);

%Technical values to build Ztilde
lag_tfp           = 8; %number of lags of TFP - cannot be zero since 1 include current TFP
lead_tfp          = 16; %number of leads of TFP
lag               = 2;  %number of lags of control variables (other structural shocks)
mpc               = 2; %max number of principal components
threshold         = -1/eps; %Remove all the NaN values

%Set Z1, Z2, and Z3 into a vector
Z                 = [Z1 Z2 Z3];

%Building Ztilde
switch 'CascaldiGarcia'
    case 'First_Principal_Component'
        for in = 1:2
            loc_start(in)   = find(Z(:,in) > threshold, 1);
        end
        loc_start       = max(loc_start) + lag;
        loc_end         = find(isnan(MUNI1Y(loc_start+1:end)),1);
        loc_end         = loc_start + loc_end - 1 - 2;
        NoZscore        = 0; %Do not standardize data before taking PC
        ZPC             = get_principal_components(Z(loc_start:loc_end-1,1:2),NoZscore);
        ZPC             = ZPC(:,1);
        ZZ              = ZPC;
    case 'CascaldiGarcia' %Get either Z1, Z2, or Z3
        iii             = 1;
        ZZ              = Z(:,iii);
        loc_start       = find(ZZ > threshold, 1) + lag;
        loc_end         = find(isnan(MUNI1Y(loc_start+1:end)),1);
        loc_end         = loc_start + loc_end - 1 - 2;
        ZZ              = ZZ(loc_start:loc_end-1);
    case 'MichiganIndex'
        ZZ              = Mich1Y;
        loc_start       = find(RESID08 > threshold, 1) + lag;
        loc_end         = find(isnan(MUNI1Y(loc_start+1:end)),1);
        loc_end         = loc_start + loc_end - 1 - 2;
        ZZ              = ZZ(loc_start:loc_end-1);
    case 'StockPrices'
        ZZ              = SP500;
        loc_start       = find(RESID08 > threshold, 1) + lag;
        loc_end         = find(isnan(MUNI1Y(loc_start+1:end)),1);
        loc_end         = loc_start + loc_end - 1 - 2;
        ZZ              = ZZ(loc_start:loc_end-1);
end

%Runniong OLS to obtain Ztilde
T                 = size(ZZ,1);
const             = ones(T,1);
X                 = const;

%Structural shocks from Ramey narrative approach, Hamilton, Romer and
%Romer, Military government spending...
controls          = [MUNI1Y, PDVMILY, HAMILTON3YP, RESID08, TAXNARRATIVE];
trend             = 1:1:length(X); %Control for the time trend
X                 = [X, trend', controls(loc_start+1:loc_end,:)];
for i = 1:lag_tfp %Add lags of TFP - When i = 1 TFP is contemporaneous
    X(:,end+1)  = TFP(loc_start+2-i:loc_end-i+1);
end
for i = 1:lead_tfp %Add leads of TFP
    X(:,end+1)  = TFP(loc_start+1+i:loc_end+i);
end
for l = 1:lag %Add lags of controls
    X           = [X controls(loc_start+1-l:loc_end-l,:) pc(loc_start+1-l:loc_end-l,1:mpc)];
end
Y                 = ZZ;
[B,zhat,Ztilde]   = quick_ols(Y,X);
for iii = 1:20
    corr(Ztilde,TFP(loc_start+1+iii:loc_end+iii));
end

%Show the graph of Ztilde - Figure(1)
Ztilde_graph = 0.05*Ztilde + .05;
figure('Position', [-1919 41 1920 963])
figure(1)
set(gcf,'color','w');
area(Time(loc_start+1:loc_end),NBERDates(loc_start+1:loc_end),...
    'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
grid on
plot(Time(loc_start+1:loc_end),Ztilde_graph,'black-','Linewidth',3)
hold off
%xlim([12 252])
ylim([.038 .060])
set(gca,'YTickLabel',[]);
lgd = legend('NBER recessions','Noise Shocks','Location',...
    'SouthOutside','Orientation','horizontal');
lgd.FontSize = 30;
legend('boxoff')

% Print figure authomatically if "export_figure1 = 1"
export_figure1 = 0;
if export_figure1 == 1
    % Create the correct path
    base_path = pwd;
    warning off
    if exist([base_path '\Figures'], 'dir')
        cd([base_path '\Figures']) %for Microsoft
    else
        cd([base_path '/Figures']) %for Mac
    end
    if exist([base_path '\Export_Fig'], 'dir')
        addpath([base_path '\Export_Fig']) %for Microsoft
    else
        addpath([base_path '/Export_Fig']) %for Mac
    end
    warning on
    export_fig(['Noise_Shocks_SPF_GDPgrowth_revisions.pdf'])
    close all
    cd(base_path) %back to the original path
end

%*************************************************************************%
%                                                                         %
%          2nd stage - Smooth Transition Local Projections                %
%                                                                         %
%*************************************************************************%

% Smooth Transition Local Projection
varlist = {'TFP','Real GDP', 'Real Consumption','Unemployment Rate','Real Wage',...
    'Hours','CPI','Real Investment','SP500','Vix','GZSpread','FFR'};
lags    = 2;
H       = 20; %irfs horizon
mpc     = 2; %max number of principal components

%standardize Ztilde to get one std dev shock
Ztilde  = Ztilde/std(Ztilde);

% Matrix of dependen variables
dep_var = [100*TFP(loc_start+1:loc_end-2) 100*RealGDP(loc_start+1:loc_end-2) ...
    100*RealCons(loc_start+1:loc_end-2)  UnempRate(loc_start+1:loc_end-2) ...
    100*RealWage(loc_start+1:loc_end-2) 100*Hours(loc_start+1:loc_end-2) ...
    100*CPIInflation(loc_start+1:loc_end-2) 100*RealInvestment(loc_start+1:loc_end-2) ...
    100*SP500(loc_start+1:loc_end-2) 100*Vix(loc_start+1:loc_end-2) ...
    100*GZSpread(loc_start+1:loc_end-2) 100*FFR(loc_start+1:loc_end-2)];

%stlp(y,x,u,fz(-1),lags,H); where y is the dep var, u is the shock, x are the controls
for kk = 1:size(dep_var,2)
    if sum(isnan(dep_var(:,kk))) > 0
        threshold = -1/eps;
        loc = find(dep_var(:,kk) > threshold, 1);
        loc_start2 = loc;
        [IR_E{kk}, IR_R{kk}, IR_L{kk}, res_uncond{kk}, Rsquared{kk}, ...
            BL{kk}, regressor{kk},SE_store{kk},SEL_store{kk},...
            tuple_store{kk}, tuple_store_conditional{kk}] = stlp(dep_var(loc_start2+1:end,kk),...
            pc(loc_start2+loc_start+1:loc_end-2,1:mpc),...
            Ztilde(loc_start2+1:end-2),ProbRecession(loc_start+loc_start2:loc_end-1-2),...
            lags,H,TFP(loc_start2+loc_start+1:loc_end-2));
    else
        [IR_E{kk}, IR_R{kk}, IR_L{kk}, res_uncond{kk}, Rsquared{kk}, ...
            BL{kk}, regressor{kk},SE_store{kk},SEL_store{kk},...
            tuple_store{kk}, tuple_store_conditional{kk}] = stlp(dep_var(:,kk),pc(loc_start+1:loc_end-2,1:mpc),...
            Ztilde(1:end-2),ProbRecession(loc_start:loc_end-1-2),...
            lags,H,TFP(loc_start+1:loc_end-2));
    end
end

nsimul = 1000;
for kk = 1:size(dep_var,2)
    tuple_depvarkk        = tuple_store{kk};
    tuple_depvarkk_cond   = tuple_store_conditional{kk};
    for hh = 1:H
        tuple_depvarkk_horizonhh        = tuple_depvarkk{hh};
        tuple_depvarkk_horizonhh_cond   = tuple_depvarkk_cond{hh};
        Y                               = tuple_depvarkk_horizonhh(:,1);
        Y_cond                          = tuple_depvarkk_horizonhh_cond(:,1);
        X                               = tuple_depvarkk_horizonhh(:,2:end);
        X_cond                          = tuple_depvarkk_horizonhh_cond(:,2:end);
        [Yboot, Xboot]                  = bb_bootstrap_LP(Y,X,nsimul,lags);
        [Yboot_cond, Xboot_cond]        = bb_bootstrap_LP(Y_cond,X_cond,nsimul,lags);
        for isimul = 1:nsimul
            B                             = Xboot(:,:,isimul)'*Xboot(:,:,isimul)\(Xboot(:,:,isimul)'*Yboot(:,isimul));
            B_cond                        = Xboot_cond(:,:,isimul)'*Xboot_cond(:,:,isimul)\(Xboot_cond(:,:,isimul)'*Yboot_cond(:,isimul));
            IRF_boot(kk,hh,isimul)        = B(1);
            IRF_boot_Exp(kk,hh,isimul)    = B_cond(1);
            IRF_boot_Rec(kk,hh,isimul)    = B_cond(1) + B_cond(2);
        end
    end
end

IRF_boot         = cumsum(IRF_boot,2);
IRF_boot         = sort(IRF_boot,3);
IRF_boot_Exp     = sort(IRF_boot_Exp,3);
IRF_boot_Rec     = sort(IRF_boot_Rec,3);
sig              = 0.05;
up_bound         = floor(nsimul*sig); % the upper percentile of bootstrapped responses for CI
low_bound        = ceil(nsimul*(1-sig)); % the lower percentile of bootstrapped responses for CI
IRF_up           = IRF_boot(:,:,up_bound);
IRF_up_Exp       = IRF_boot_Exp(:,:,up_bound);
IRF_up_Rec       = IRF_boot_Rec(:,:,up_bound);
IRF_low          = IRF_boot(:,:,low_bound);
IRF_low_Exp      = IRF_boot_Exp(:,:,low_bound);
IRF_low_Rec      = IRF_boot_Rec(:,:,low_bound);

% Build a table for the Rsquared
% This R-squared has to be interpreted as the variance explained by noise
% shocks of macroeconomic variables at each specific horizon
for kkk = 1:size(dep_var,2) %Raws are time horizons, Columns are variables.
    Rsquared_Table(:,kkk) = Rsquared{kkk}';
end
varlist;
Rsquared_Table;

%Impulse Response Functions using Local Projection
nvar     = length(varlist);
n_row    = 2;
n_col    = ceil(nvar/n_row); %plus one for Vix
figure('Position',[1 41 1920 963])
set(gcf,'color','w');

for j = 1:length(varlist)
    s = subplot(n_row,n_col,j);
    hold on
    if j >= 100 %Be careful if j >= 1 No variables are cumulated!
        %             h1 = plot([1:H]',IR_L{j}+1.64*SEL_store{j}, '--k','linewidth', 2);
        %             h2 = plot([1:H]',IR_L{j}-1.64*SEL_store{j}, '--k','linewidth', 2);
        h1  = plot([1:H]',IRF_low(j,:), '--k','linewidth', 2);
        h1  = plot([1:H]',IRF_up(j,:), '--k','linewidth', 2);
        q   = plot([1:H]',IR_L{j}, 'k', 'linewidth', 3);
    else
        % h1 = plot([1:H]',cumsum(IR_L{j}+1.96*SEL_store{j}), '.k','linewidth', 2);
        % h2 = plot([1:H]',cumsum(IR_L{j}-1.96*SEL_store{j}), '.k','linewidth', 2);
                     h1  = plot([1:H]',IRF_low(j,:), '--k','linewidth', 2);
                     h1  = plot([1:H]',IRF_up(j,:), '--k','linewidth', 2);
                     q   = plot([1:H]',cumsum(IR_L{j}), '-k', 'linewidth', 3);
        %             h = plot([1:H]',cumsum(IR_R{j}), '--b','linewidth', 3);
        %l = plot([1:H]',cumsum(IR_L{j}), '-ok','linewidth', 3);
    end
        plot([1:H]', 0*[1:H]', ':k');
        set(gca,'TickLabelInterpreter','latex')
        title(varlist{j},'interpreter', 'latex', 'fontsize', 14);
    if j == 1
        xlabel('Quarter','interpreter','latex','fontsize',12);
        ylabel('\% deviation from s.s.','interpreter','latex','fontsize',12);
    end
    set(s, 'xlim', [1,H], 'ylim', ylim );
end

%l = legend([q h l],{'Expansion','Recession','Linear'},'interpreter','latex');
%set(l, 'box','off', 'FontSize',30,'Orientation','horizontal','Position',[0.3 0.015 0.40 0.01]);

% Print figure automatically if "export_figure1 = 1"
export_figure2 = 0;
if export_figure2 == 1
    % Create the correct path
    base_path = pwd;
    warning off
    if exist([base_path '\Figures'], 'dir')
        cd([base_path '\Figures']) %for Microsoft
    else
        cd([base_path '/Figures']) %for Mac
    end
    if exist([base_path '\Export_Fig'], 'dir')
        addpath([base_path '\Export_Fig']) %for Microsoft
    else
        addpath([base_path '/Export_Fig']) %for Mac
    end
    warning on
    export_fig(['STLP_IRFs_Noise_Shocks_SPF_GDPgrowth_revisions_logdiff_bootstrap.pdf'])
    close all
    cd(base_path) %back to the original path
end














