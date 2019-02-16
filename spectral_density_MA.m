function [sdensity, sdensity_up, sdensity_low, ...
      sdensity_up2, sdensity_low2, sdensity_median, sdensity_boot, period] ...
      = spectral_density_MA(IRF,IRF_boot,sig,sig2,showfig)

% This function takes IRF and IRF_boot as truncated approximation of MA
% representation of endogenous variables conditional on a shock. Given
% MA(H), it generates spectral density for both point estimate and
% bootstrap estimations. sig and sig2 determines the confidence intervals.
% showfig = 1, plot spectral density. otherwide do not plot it. 

% Evaluate Spectral Density
[sdensity_boot]                 = spectrum(IRF_boot); % bootstrap 
[sdensity, period(:,1)]         = spectrum(IRF); %point estimate

% Normalization
sdensity_boot                   = sdensity_boot.*1/(2*sum(sdensity));
sdensity                        = sdensity.*1/(2*sum(sdensity));

% Determine Confidence intervals with sig and sig2
sdensity_up(:,1)                = quantile(sdensity_boot',1-sig);
sdensity_low(:,1)               = quantile(sdensity_boot',sig);
sdensity_up2(:,1)               = quantile(sdensity_boot',1-sig2);
sdensity_low2(:,1)              = quantile(sdensity_boot',sig2);
sdensity_median(:,1)            = quantile(sdensity_boot',.5);

% Plot figure if showfig = 1
per_start = 10;
per_end   = 200;
if showfig == 1
      %plot spectral density of AR(1) from LP
      hfig        = findobj('type','figure');
      nfig        = length(hfig);
      figure(1+nfig)
      set(gcf,'Position',[-1919 41 1920 963])
      set(gcf,'color','w');
      plot(period(per_start:per_end)',sdensity(per_start:per_end),'-r','LineWidth',2); hold on; %step dependent
      if IRF_boot == 0
            lgd = legend('Point Estimate','Location',...
                  'South','Orientation','horizontal');
      else
            plot(period(per_start:per_end)',sdensity_median(per_start:per_end),'-k','LineWidth',3); hold on;
            plot(period(per_start:per_end)',sdensity_up(per_start:per_end),'--k','LineWidth',2); hold on;
            plot(period(per_start:per_end)',sdensity_low(per_start:per_end),'--k','LineWidth',2); hold on; %the point estimate is not included in the CI, is it because we don't correct for the bias in the LP?
            lgd = legend('Point Estimate','Median','Location',...
                  'South','Orientation','horizontal');
      end
      lgd.FontSize = 24;
      legend('boxoff')
      title('Spectral Density','fontsize',26)
      axis tight
end





end