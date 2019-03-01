fprintf('\n')
fprintf('\n')
disp('-------------------------------------------------------------------')
disp('                      Technical Parameters LP                      ')
disp('-------------------------------------------------------------------')
for iw = 1: length(which_Z)
if strcmp(which_Z{iw},'1') == 1
      disp('Survey                        SPF Real GDP Growth')
elseif strcmp(which_Z{iw},'2') == 1
      disp('Survey                        SPF Nominal GDP Growth')
elseif strcmp(which_Z{iw},'3') == 1
      disp('Survey                        SPF Real Consumption Growth')
elseif strcmp(which_Z{iw},'4') == 1
      disp('Survey                        SPF Industrial Production Growth')
elseif strcmp(which_Z{iw},'5') == 1
      disp('Survey                        SPF Real Total Investment Growth')
elseif strcmp(which_Z{iw},'6') == 1
      disp('Survey                        SPF Consumer Price Index')
elseif strcmp(which_Z{iw},'7') == 1
      disp('Survey                        Michigan Index Confidence')
end
end
disp(['Number of lags 1st Step       ',num2str(lags)])
disp(['Number of leads               ',num2str(leads)])
disp(['Number of lags in LP          ',num2str(lags_LP)])
disp(['Trend                         ',which_trend])
disp(['IRFs horizon                  ',num2str(H)])
disp(['Number of PCs 1st Step        ',num2str(nPC_first)])
disp(['Number of PCs in LP           ',num2str(nPC_LP)])
disp(['Number of simulations         ',num2str(nsimul)])
disp(['Number of regressors in LP    ',num2str(floor(mean(nREGkk(1,:,1))))])
if diff_LP == 1 && control_pop == 1
      disp('Other                         Dependent Variables are differentiated and per capita')
elseif diff_LP == 1 && control_pop == 0
      disp('Other                         Dependent Variables are differentiated')
elseif diff_LP == 0 && control_pop == 1
      disp('Other                         Dependent Variables are per capita')
else
      disp('Other                         Dependent Variables are in levels')
end
disp('-------------------------------------------------------------------')
fprintf('\n')
for ishock = 1:length(which_shock)
      disp([which_shock{ishock}])
      for ii = 1:length(varlist)
            disp(['  - Degrees of freedom for ',varlist{ii},' are on average ',num2str(floor(mean(DFkk(ii,:,ishock)))),' and first observation is in ',num2str(Time(loc_start(ii,ishock)))])
      end
      fprintf('\n')
end
