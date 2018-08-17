function plot_IRF_lp_unconditional(varlist,IRF_E_low,IRF_E_up,IRF_E,...
IRF_R_low,IRF_R_up,IRF_R,H,plot2,n_row,unique)

if plot2 == 1
      %Impulse Response Functions using Local Projection - Figure
      if unique == 1
            nvar     = length(varlist);
            n_col    = ceil(nvar/n_row); %plus one for Vix
            figure('Position',[1 41 1920 963])
            set(gcf,'color','w');
      end      
      for j = 1:length(varlist)
            if unique == 1
                  s = subplot(n_row,n_col,j);
            else
                  figure(j)
                  nvar     = length(varlist);
                  n_col    = ceil(nvar/n_row); %plus one for Vix
                  figure('Position',[1 41 1920 963])
                  set(gcf,'color','w');
            end
            hold on
            plot([0:H-1]',IRF_E_low(j,:), '--b','linewidth', 1);
            plot([0:H-1]',IRF_E_up(j,:), '--b','linewidth', 1);
            plot([0:H-1]',IRF_E(j,:), '-b', 'linewidth', 3);
            plot([0:H-1]',IRF_R_low(j,:), '--r','linewidth', 1);
            plot([0:H-1]',IRF_R_up(j,:), '--r','linewidth', 1);
            plot([0:H-1]',IRF_R(j,:), '-r', 'linewidth', 3);
            plot([0:H-1]',0*[1:H]',':k');
            set(gca,'TickLabelInterpreter','latex')
            title(varlist{j},'interpreter', 'latex', 'fontsize', 14);
            if unique == 1 && j == 1
                  xlabel('Quarter','interpreter','latex','fontsize',12);
                  ylabel('\% deviation from s.s.','interpreter','latex','fontsize',12);
            elseif unique == 0
                  xlabel('Quarter','interpreter','latex','fontsize',12);
                  ylabel('\% deviation from s.s.','interpreter','latex','fontsize',12);
            end
            if unique == 1
                  set(s,'xlim',[1,H],'ylim', ylim);
            end
            
      end
end