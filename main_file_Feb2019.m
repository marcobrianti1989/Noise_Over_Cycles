%*************************************************************************%
%    Main
%
%    NOTE describe variables (especially SHOCKS) in dataset
%
%    last change 02/14/2019
%
%    Code by Brianti, Marco e Cormun, Vito
%*************************************************************************%
clc
clear
close all

% Technical Parameters
lags                = 4;             % Number of lags in the first step (deriving Ztilde)
leads               = 16;            % Number of leads in the first step (deriving Ztilde)
H                   = 20;            % IRFs horizon
lags_LP             = 2;             % Number of lags in the Local Projection
which_trend         = 'quadratic' ;  % BPfilter, HPfilter, linear, quadratic for Local Projection
which_Z             = '1';           % Which Forecast Revision: RGDP, NGDP, RCONS, INDPROD, RINV
which_shock         = {'Sentiment','Tech'}; % Tech, News
diff_LP             = 0;             % LP in levels or differences
nPC                 = 3;             % Number of Principal Components
norm_SHOCK          = 1;             % Divide shock over its own variance
printIRFs           = 1;             % Print IRFs
printVD             = 1;             % Print Variance Decompositions

% Read main dataset
filename                    = 'main_file';
sheet                       = 'Sheet1';
range                       = 'B1:CX300';
do_truncation               = 0; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data2(filename, sheet, range, do_truncation);
dataset                     = [dataset; NaN(leads,size(dataset,2))]; % Adding some NaN at the end for technical purposes

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
create_Z

% Define Variables
[TFP_trunc, trunc1, trunc2] = truncate_data(TFP);
TFPBP                       = bpass(TFP_trunc,4,32);
TFPBP                       = [TFPBP; NaN(length(TFP) - length(TFPBP),1)];
dTFP                        = [NaN; diff(TFP)];
PC                          = [PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9];
PC                          = PC(:,1:nPC);
SHOCKS_NARRATIVE            = [MUNI1Y PDVMILY HAMILTON3YP RESID08 TAXNARRATIVE];

% Defibe dependent variable
eval(['Z = Z', which_Z,';']);
Y                           = Z;

% Define Regressors and Dependent Variable
X_contemporaneous           = [TFPBP SHOCKS_NARRATIVE];
X_lag                       = [TFPBP PC SHOCKS_NARRATIVE];
X_lead                      = TFPBP;

% Control Regression
[~, Zhat, Ztilde, regressor] = lead_lag_matrix_regression(Y,X_lead,...
      leads,X_lag,lags,X_contemporaneous);

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
control_pop = 0; % Divide GDP, Cons, Hours, Investment over population
if control_pop == 1
      RealGDP                 = RealGDP - Population;
      RealCons                = RealCons - Population;
      RealInvestment          = RealInvestment - Population;
      Hours                   = Hours + Employment - Population;
      RealInventories         = RealInventories - Population;
      RealSales               = RealSales - Population;
      SP500                   = SP500 - GDPDefl;
      Spread                  = MoodySpreadBaa - MoodySpreadAaa;
else
      SP500                   = SP500 - GDPDefl;
      Spread                  = MoodySpreadBaa - MoodySpreadAaa;
end

% Define Dependent Variables
varlist          = {'RealGDP','RealCons','RealInvestment','RealInventories',...
      'Credit2GDP','Spread'};

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
PC               = PC(1+lags:end-leads,:);

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
      IRFname = ['IRFs_',char(which_shock(is)),'Shock_lags',num2str(lags),'_','_leads',num2str(leads),'_lagsLP',num2str(lags_LP),'_trend',which_trend,'_Z',num2str(which_Z),'_diff',num2str(diff_LP),'_nPC',num2str(nPC),'.pdf'];
      VDname  = ['VD_',char(which_shock(is)),'Shock_lags',num2str(lags),'_','_leads',num2str(leads),'_lagsLP',num2str(lags_LP),'_trend',which_trend,'_Z',num2str(which_Z),'_diff',num2str(diff_LP),'_nPC',num2str(nPC)','.pdf'];
      
      % Normilize Variance of SHOCK
      if norm_SHOCK == 1
            sdSHOCK          = nanstd(SHOCK);
            SHOCK            = SHOCK/sdSHOCK;
      end
      
      %Initializating the loop
      for kk = 1:size(dep_var,2)
            % Define inputs for local_projection
            varnamekk                   = varlist{kk};
            disp(['Projecting ',varnamekk])
            fprintf('\n')
            depvarkk                    = dep_var(:,kk);
            [~, loc_start, loc_end]     = truncate_data([depvarkk SHOCK PC]);
            depvarkk                    = depvarkk(loc_start:loc_end);
            SHOCKkk                     = SHOCK(loc_start:loc_end);
            pckk                        = PC(loc_start:loc_end,:);
            % Run local_projection
            [IR{kk},res{kk},tuple{kk},VD{kk}] = ...
                  local_projection(depvarkk,pckk,SHOCKkk,lags_LP,H,which_trend);
            if diff_LP == 0
                  IRF(kk,:) = IR{kk};
            else
                  IRF(kk,:) = cumsum(IR{kk});
            end
            VDkk(kk,:) = VD{kk};
            % Initiate bootstrap
            nsimul         = 1000;
            tuplekk        = tuple{kk};
            for hh = 1:H
                  tuplekkhh = tuplekk{hh}; % Fix a specific horizon
                  Y                             = tuplekkhh(:,1);
                  X                             = tuplekkhh(:,2:end);
                  [Yboot, Xboot]                = bb_bootstrap_LP(Y,X,nsimul,lags_LP);
                  for isimul = 1:nsimul
                        B                       = Xboot(:,:,isimul)'*Xboot(:,:,isimul)\...
                              (Xboot(:,:,isimul)'*Yboot(:,isimul));
                        IRF_boot(kk,hh,isimul)  = B(1);
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
            IRF_up(j,:)   = quantile(squeeze(IRF_boot(j,:,:))',1-sig);
            IRF_up2(j,:)  = quantile(squeeze(IRF_boot(j,:,:))',1-sig2);
            IRF_low(j,:)  = quantile(squeeze(IRF_boot(j,:,:))',sig);
            IRF_low2(j,:) = quantile(squeeze(IRF_boot(j,:,:))',sig2);
      end
      
      % Plot IRFs
      plot_IRF(varlist,IRF_low,IRF_low2,IRF_up,IRF_up2,IRF,H,printIRFs,IRFname); %change this function
      % Plot VD
      plot_IRF(varlist,VDkk,VDkk,VDkk,VDkk,VDkk,H,printVD,VDname)
end

% Technical Parameters
tech_info_table;
asd

%*************************************************************************%
%                                                                         %
%                        Test Cyclicality - Canova                        %
%                                                                         %
%*************************************************************************%

% Generate IRF on AR(1) from LP
rho    = 0;    % rho = 0, Hnull: flat spectral density, otherwise rho should be estimated from data as im BG
T      = 1000; % Asyntotic (Maybe should be the same length) !!!
Yar(1) = 0;
for t = 2:T
      Yar(t) = rho*Yar(t-1) + randn;
end
% Recover shocks
[coef,~,Zar] = regress(Yar(2:end)',Yar(1:end-1)');
% Run LP
[IRFar,res,tuplear] = local_projection(Yar(2:end)',zeros(T-1,1),Zar,0,H,'none');
for hh = 1:H
      tuplearhh = tuplear{hh}; % Fix a specific horizon
      Y                             = tuplearhh(:,1);
      X                             = tuplearhh(:,2:end);
      [Yboot, Xboot]                = bb_bootstrap_LP(Y,X,nsimul,0);
      for isimul = 1:nsimul
            B                       = Xboot(:,:,isimul)'*Xboot(:,:,isimul)\...
                  (Xboot(:,:,isimul)'*Yboot(:,isimul));
            IRFar_boot(1,hh,isimul)  = B(1);
      end
end
% Plot LP under AR(1)
IRFar_up    = quantile(squeeze(IRFar_boot(1,:,:))',1-sig);
IRFar_up2   = quantile(squeeze(IRFar_boot(1,:,:))',1-sig2);
IRFar_low   = quantile(squeeze(IRFar_boot(1,:,:))',sig);
IRFar_low2  = quantile(squeeze(IRFar_boot(1,:,:))',sig2);
hfig        =  findobj('type','figure');
nfig        = length(hfig);
figure(nfig); %plot spectral density and its CI against the AR(1) counterpart
print       = NaN;
nameAR      = {'AR1'};
plot_IRF(nameAR,IRFar_low,IRFar_low2,IRFar_up,IRFar_up2,IRFar,H,print,name); %change this function

[sdensityar] = spectrum(IRFar_boot);
[sdensityar_pe, period] = spectrum(IRFar); %point estimate

%plot spectral density of AR(1) from LP
figure(4+cc);
sdensityar_up   = quantile(sdensityar',1-sig);
sdensityar_low  = quantile(sdensityar',sig);
sdensityar_ave  = quantile(sdensityar',.5);
plot(period(10:200)',sdensityar_pe(10:200),'-r','LineWidth',2); hold on; %step dependent
plot(period(10:200)',sdensityar_ave(10:200),'-b','LineWidth',3); hold on;
plot(period(10:200)',sdensityar_up(10:200),'--b','LineWidth',2); hold on;
plot(period(10:200)',sdensityar_low(10:200),'--b','LineWidth',2); hold on; %the point estimate is not included in the CI, is it because we don't correct for the bias in the LP?

% Construct spectrum from data
[sdensity] = spectrum(IRF_boot(1,:,:));
[sdensity_pe(:,cc), period] = spectrum(IRF(1,:)); %point estimate
sdensity_up(cc,:)   = quantile(sdensity',1-sig);
sdensity_low(cc,:)  = quantile(sdensity',sig);
sdensity_ave(cc,:)  = quantile(sdensity',.5);
%Normalize
sdensity_up(cc,:) = sdensity_up(cc,:).* 1/(2*sum(sdensity_up(cc,:)));
sdensity_low(cc,:) = sdensity_low(cc,:).* 1/(2*sum(sdensity_low(cc,:)));
sdensity_ave(cc,:) = sdensity_ave(cc,:).* 1/(2*sum(sdensity_ave(cc,:)));
%Compute average spectral density, D1, around the peak  and average
%spectral density around the trough, D2
lpeak_lower   = 24; %should be adjusted with steps and IRF horizon
lpeak_upper   = 26;
ltrough_lower = 58;
ltrough_upper = 60;
D1 = mean(sdensity(find(period>lpeak_upper,1,'last'):find(period>lpeak_lower,1,'last'),:),1);
D2 = mean(sdensity(find(period>ltrough_upper,1,'last'):find(period>ltrough_lower,1,'last'),:),1);
D  = D1./D2;

D1ar = mean(sdensityar(find(period>lpeak_upper,1,'last'):find(period>lpeak_lower,1,'last'),:));
D2ar = mean(sdensityar(find(period>ltrough_upper,1,'last'):find(period>ltrough_lower,1,'last'),:));
Dar = D1ar./D2ar;

Diff_D = D - Dar;
pval(cc) = 1 - length(find(Diff_D>0))/nsimul %results seems to favor white noise against ar(1)
%end end for cc = 1:2
var_list = {'$S_{GDP}$ to Sentiment','$S_{GDP}$ to Technology'};
for cc = 1:2
      %plot spectral density and its CI against the AR(1) counterpart
      figure(500);
      a = subplot(1,2,cc);
      p = plot(period(10:200)',sdensity_pe(10:200,cc),'-r','LineWidth',2); hold on; %step dependent
      m = plot(period(10:200)',sdensity_ave(cc,10:200),'-b','LineWidth',3); hold on;
      plot(period(10:200)',sdensity_up(cc,10:200),'--b','LineWidth',2); hold on;
      plot(period(10:200)',sdensity_low(cc,10:200),'--b','LineWidth',2); hold on; %the point estimate is not included in the CI, is it because we don't correct for the bias in the LP?
      t = plot(period(10:200)',sdensityar_pe(10:200),'k','LineWidth',3);
      title(var_list{cc},'interpreter', 'latex', 'fontsize', 16);
      if cc == 1
            xlabel('Periodicity','fontsize',12);
      end
      if cc == 2
            l=legend([p m t],{'Point estimate','Median','Null'},'Location', 'NorthEast','interpreter','latex');
            set(l, 'box','off', 'FontSize',14,'Orientation','horizontal','Position',[0.35179282868526 -0.00717088619666764 0.400000000000001 0.0617721521401707]);
      end
end




%% Trash
% %% Multivariate test
% % it should embed the idea that there should be a pick in both spectral
% % densities of the series and their coherence should be high at the
% % frequency considered
%
% %Spectral PCA
% % standardize IRFs
% IRF_boot_z = IRF_boot./sum(IRF_boot.^2,2);
%
% % compute spectrum
% [spect, periodicity] = spectrum(IRF_boot_z);
%
% % compute PCA
% % take 1st pc and construct D
%
% [V,lamb] = eig(sigma); % diag(lamb) are the eigenvalues, V the eigenvectors
% V        = real(V);
% lamb     = real(lamb);
% pc       = data*V./n;
%
% % %Theoretical AR1-IRF
% % IRF_ar1 = 1;
% % for j = 2: H
% %       IRF_ar1(j,:) = rho^j;
% % end
% % [sdensity_ar1] = spectrum(IRF_ar1); %I normalize s.t. the spectral density evaluated btwn 0 and pi is equal to .5 - CHECK
% % plot(sdensity_ar1)
% %
% % figure(5); %plot spectral density and its CI against the AR(1) counterpart
% % sdensity_up   = quantile(sdensity',1-sig);
% % sdensity_low  = quantile(sdensity',sig);
% % sdensity_ave  = quantile(sdensity',.5);
% % plot(period(10:200)',sdensity_pe(10:200),'-r','LineWidth',2); hold on; %step dependent
% % plot(period(10:200)',sdensity_ave(10:200),'-b','LineWidth',3); hold on;
% % plot(period(10:200)',sdensity_up(10:200),'--b','LineWidth',2); hold on;
% % plot(period(10:200)',sdensity_low(10:200),'--b','LineWidth',2); hold on; %the point estimate is not included in the CI, is it because we don't correct for the bias in the LP?
% % plot(period(10:200)',sdensityar_pe(10:200),'k','LineWidth',1.25);
% % xlabel('Frequency','fontsize',20);
% % grid on
% %
% % %Compute average spectral density, D1, around the peak  and average
% % %spectral density around the trough, D2
% % lpeak_lower   = 25; %should be adjusted with steps and IRF horizon
% % lpeak_upper   = 35;
% % ltrough_lower = 40;
% % ltrough_upper = 50;
% % D1 = mean(sdensity(find(period>lpeak_upper,1,'last'):find(period>lpeak_lower,1,'last'),:),1);
% % D2 = mean(sdensity(find(period>ltrough_upper,1,'last'):find(period>ltrough_lower,1,'last'),:),1);
% % D  = D1./D2;
% % D1_ar1 = mean(sdensity_ar1(find(period>lpeak_upper,1,'last'):find(period>lpeak_lower,1,'last'),:));
% % D2_ar1 = mean(sdensity_ar1(find(period>ltrough_upper,1,'last'):find(period>ltrough_lower,1,'last'),:));
% % D_ar1 = D1_ar1./D2_ar1;
% % pval = length(find(D<=D_ar1))/nsimul;
% % disp(['P-value univariate:    ',num2str(pval)])
% %
% % %Need to perform LP on the AR(1)









