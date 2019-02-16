function [IRFar, IRFar_boot] = generate_AR1_LP_IRF(rho,sigm,T,H,nsimul,sig,sig2,showfig)

% This function simulate T observations of an AR1 with persistence rho and
% variance sigm.
% Once simulated the series, it generates MA(H) representation via Local
% Projection and evaluate confidence interval via bootstrap with nsimul
% simulations. sig and sig2 are confidence intervals parameters.
% Showfig = 1 plot the figures, otherwise it does not plot.

% Generate IRF on AR(1) from LP
Yar(1) = 0;
for t = 2:T
      Yar(t) = rho*Yar(t-1) + sigm*randn;
end
% Recover shocks for LP
[~,~,Zar] = regress(Yar(2:end)',Yar(1:end-1)');
% Run LP
[IRFar,~,tuplear] = local_projection(Yar(2:end)',zeros(T-1,1),Zar,0,H,'none');
% Inference via Bootstrap
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
% Use sig and sig2 to identify confidence intervals
IRFar_up    = quantile(squeeze(IRFar_boot(1,:,:))',1-sig);
IRFar_up2   = quantile(squeeze(IRFar_boot(1,:,:))',1-sig2);
IRFar_low   = quantile(squeeze(IRFar_boot(1,:,:))',sig);
IRFar_low2  = quantile(squeeze(IRFar_boot(1,:,:))',sig2);

% Plot (if showfig = 1) LP under AR(1)
if showfig == 1
      hfig        = findobj('type','figure');
      nfig        = length(hfig);
      figure(1+nfig)
      set(gcf,'Position',[1 41 1920 963])
      set(gcf,'color','w');
      nameAR      = {['IRF of AR(1) with rho = ',num2str(rho)]};
      name        = ['IRF of AR(1) with \rho = ',num2str(rho)];
      hold on
      plot([0:H-1]',IRFar_low, '--k','linewidth', 1);
      plot([0:H-1]',IRFar_up, '--k','linewidth', 1);
      plot([0:H-1]',IRFar_low2, '--k','linewidth', 2);
      plot([0:H-1]',IRFar_up2, '--k','linewidth', 2);
      plot([0:H-1]',IRFar, '-k', 'linewidth', 3);
      plot([0:H-1]',0*[1:H]',':k');
      set(gca,'TickLabelInterpreter','latex')
      xlabel('Quarter','interpreter','latex','fontsize',20);
      ylabel('\% deviation from s.s.','interpreter','latex','fontsize',18);
      title(['IRF AR(1) with persistence ',num2str(rho)],'fontsize',26)
      axis tight
end