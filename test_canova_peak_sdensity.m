function pval = test_canova_peak_sdensity(lpeak_lower,lpeak_upper,ltrough_lower,...
      ltrough_upper,sdensity,sdensityar,period,nsimul)

% Test Canova (1996) also used by Beaudry, Galizia, Portier (2019, AER)
% Given a spectral density, set D1 equal to the average SD explained over
% business cycle frequencies (lpeak_lower,lpeak_upper). Set D2 equal to the average SD explained over
% the long run (ltrough_lower,ltrough_upper). Take the ratio D = D1/D2.
% Compare D between empirical IRFs and simulated AR(1) process: 
% Diff_D = D - Dar

% Empirical D
D1 = mean(sdensity(find(period>lpeak_upper,1,'last'):find(period>lpeak_lower,1,'last'),:),1);
D2 = mean(sdensity(find(period>ltrough_upper,1,'last'):find(period>ltrough_lower,1,'last'),:),1);
D  = D1./D2;

% Simulated AR(1) D
D1ar = mean(sdensityar(find(period>lpeak_upper,1,'last'):find(period>lpeak_lower,1,'last'),:));
D2ar = mean(sdensityar(find(period>ltrough_upper,1,'last'):find(period>ltrough_lower,1,'last'),:));
Dar = D1ar./D2ar;

% Test. H0: Empirical IRFs are AR(1). If pvalue is close to zero. Reject
% H0. Empirical IRFs are cyclical. 
Diff_D = D - Dar;
pval   = 1 - length(find(Diff_D>0))/nsimul; 


end


