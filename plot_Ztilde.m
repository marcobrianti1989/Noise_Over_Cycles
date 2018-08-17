function plot_Ztilde(Ztilde,Time,NBERDates,loc_start,loc_end,plotYes)

if plotYes == 1
      Ztilde_graph   = 0.8*Ztilde + .05;
      Time           = Time(loc_start:loc_end);
      NBERDates      = NBERDates(loc_start:loc_end);
      figure('Position', [-1919 41 1920 963])
      figure(1)
      set(gcf,'color','w');
      area(Time,NBERDates,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
      hold on
      grid on
      plot(Time,Ztilde_graph,'black-','Linewidth',3)
      hold off
      %xlim([12 252])
      ylim([.038 .060])
      set(gca,'YTickLabel',[]);
      lgd = legend('NBER recessions','Noise Shocks','Location',...
            'SouthOutside','Orientation','horizontal');
      lgd.FontSize = 30;
      legend('boxoff')
end

end
