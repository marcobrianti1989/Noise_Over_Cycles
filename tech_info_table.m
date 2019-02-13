fprintf('\n')
fprintf('\n')
fprintf('\n')
fprintf('\n')
disp('------------------------------------------------')
disp('             Technical Parameters               ')
disp('------------------------------------------------')
fprintf('\n')
if strcmp(which_Z,'1') == 1
      disp('Survey                  SPF Real GDP Growth')
elseif strcmp(which_Z,'2') == 1
      disp('Survey                  SPF Nominal GDP Growth')
elseif strcmp(which_Z,'3') == 1
      disp('Survey                  SPF Real Consumption Growth')
elseif strcmp(which_Z,'4') == 1
      disp('Survey                  SPF Industrial Production Growth')
elseif strcmp(which_Z,'5') == 1
      disp('Survey                  SPF Real Total Investment Growth')
elseif strcmp(which_Z,'6') == 1
      disp('Survey                  SPF Consumer Price Index')
elseif strcmp(which_Z,'7') == 1
      disp('Survey                  Michigan Index Confidence')
end
disp(['Number of lags          ',num2str(lags)])
disp(['Number of leads         ',num2str(leads)])
disp(['Number of lags in LP    ',num2str(lags_LP)])
disp(['Trend                   ',which_trend])
disp(['IRFs horizon            ',num2str(H)])
disp(['Number of PCs           ',num2str(nPC)])
if diff_LP == 1
      disp('Other                   Dependent Variables are differentiated')
else
      disp('Other                   Dependent Variables are in levels')
end
fprintf('\n')
fprintf('\n')
fprintf('\n')
fprintf('\n')