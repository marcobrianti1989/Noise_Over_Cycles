fprintf('\n')
fprintf('\n')
disp('-------------------------------------------------------------------')
disp('                Technical Parameters Recursive VAR                 ')
disp('-------------------------------------------------------------------')
disp(['Number of lags in VAR                   ',num2str(nlags)])
disp(['Trend                                   ',which_trend])
disp(['Specific Bootstrap Technique            ',which_boot])
disp(['IRFs horizon                            ',num2str(H)])
disp(['Number of simulations                   ',num2str(nsimul)])
disp(['Confidence intervals are                ',num2str((1-sig1*2)*100),'%-',num2str((1-sig2*2)*100),'%'])
if control_pop == 1
      disp('Other                                   Dependent Variables are per capita')
end
for i = 1:length(which_shocks)
      disp([shocknames{i},' Position is             ',num2str(which_shocks(i))])
end
disp('-------------------------------------------------------------------')


