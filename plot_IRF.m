function plot_IRF(varlist,IRF_low,IRF_low2,IRF_up,IRF_up2,IRF,H,...
      exp_fig,fig_name)

% Technical parameters, you may want to change depending on the number of Nvar
if length(varlist) == 1
      unique     = 0; % if unique = 1 plot IRFs together, if = 1 plot each IRF separately
      n_row      = NaN;
elseif length(varlist) <= 3
      unique     = 1;
      n_row      = 1; % number of rows in the figure
elseif length(varlist) <= 6
      unique     = 1;
      n_row      = 2; % number of rows in the figure
elseif length(varlist) > 6
      unique     = 1;
      n_row      = 3; % number of rows in the figure
end

% IRFs are now percentage deviations
IRF_low    = IRF_low*100;
IRF_low2   = IRF_low2*100;
IRF_up     = IRF_up*100;
IRF_up2    = IRF_up2*100;
IRF        = IRF*100;

% Periods for shadowed areas
periods = 0:H-1;

%Impulse Response Functions using Local Projection - Figure
hfig        = findobj('type','figure');
nfig        = length(hfig);
if unique == 1
      nvar     = length(varlist);
      n_col    = ceil(nvar/n_row); %plus one for Vix
end
for j = 1:length(varlist)
      if unique == 1
            figure(nfig+1)
            s = subplot(n_row,n_col,j);
      else
            figure(nfig+j)
            nvar     = length(varlist);
            n_col    = ceil(nvar/n_row); %plus one for Vix
      end
      set(gcf,'Position',[1 41 1920 963])
      set(gcf,'color','w');
      hold on
      
        x2 = [periods, fliplr(periods)];
        % The smaller CI
        inBetween = [IRF_low(j,:),fliplr(IRF_up(j,:))];
        hh2 = fill(x2, inBetween, [0.55 0.55 0.55],'LineStyle','none');
        set(hh2,'facealpha',.5)
        % The bigger CI
        inBetween = [IRF_low2(j,:),fliplr(IRF_up2(j,:))];
        hh1 = fill(x2, inBetween, [0.15 0.15 0.15],'LineStyle','none');
        set(hh1,'facealpha',.5)        

        
%       plot([0:H-1]',IRF_low(j,:), '--k','linewidth', 1);
%       plot([0:H-1]',IRF_up(j,:), '--k','linewidth', 1);
%       plot([0:H-1]',IRF_low2(j,:), '--k','linewidth', 2);
%       plot([0:H-1]',IRF_up2(j,:), '--k','linewidth', 2);
        plot([0:H-1]',IRF(j,:), '--k', 'linewidth', 2,'color','r');
        plot([0:H-1]',0*[1:H]','-k','color','b');


      set(gca,'TickLabelInterpreter','latex')
      title(varlist{j},'interpreter', 'latex', 'fontsize', 26);
      if unique == 1 && j == 1
            xlabel('Quarter','interpreter','latex','fontsize',20);
            ylabel('\% deviation from s.s.','interpreter','latex','fontsize',18);
      elseif unique == 0
            xlabel('Quarter','interpreter','latex','fontsize',20);
            ylabel('\% deviation from s.s.','interpreter','latex','fontsize',18);
      end
      axis tight
      
      % Print Fig
      if exp_fig == 1
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
            export_fig(fig_name)
            cd(base_path) %back to the original path
      end
end
