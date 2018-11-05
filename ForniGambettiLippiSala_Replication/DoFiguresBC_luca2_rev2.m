

% if kilian_flag == 0;
%     irf_chol = CI_chol(varseries,:,:,:);
%     irfsb = CI(varseries,:,:,:);
% else
%     irf_chol = CI_chol;%(varseries,:,:,:);%se usi Kilian
%     irfsb = CI;% %se usi Kilian: 
% end

%irf_plot = irf(varseries,:,:);

ff = ff+1;
figure(ff);
lo=0;
tt=[ind1  ind];
for t= [ind1 indControl ind]
    for j=1:2,
        temp = squeeze(irf_chol(t,tt(j),1:ll,:));
        TEMP = irf_plot(t,tt(j),1:ll);
        TEMP2 = irf_plot_rev(t,tt(j),1:ll);
        lo=lo+1;
        subplot(3,2,lo);
        te2 = prctile(temp',[5])';
        te3 = prctile(temp',[95])';
        te4 = prctile(temp',[16])';
        te5 = prctile(temp',[84])';
        te_med = prctile(temp',[50])';
        b=fill([0:ll-1 ll-1:-1:0],[te2' flipud(te3)'],col1);set(b,'EdgeAlpha',0);hold on
        b=fill([0:ll-1 ll-1:-1:0],[te4' flipud(te5)'],col2);set(b,'EdgeAlpha',0);
        a = plot(0:ll-1,[squeeze(TEMP) squeeze(TEMP2) 0*te_med ]);
        set(a(1),'Color','k','LineWidth',2);
        set(a(2),'Color','r','LineWidth',2,'LineStyle','--');
        set(a(end),'Color','k');
        TE = [te2 te3 te4 te5];
        TE2(:,:,:,j) = TE;
        %c = plot(0:ll-1,TE,'k');
        if lo==1
            c=title('Surprise');set(c,'FontSize',12);
            ylabel(ti{t},'FontSize',12);
        elseif lo==2
            c=title('Signal');set(c,'FontSize',12);
        elseif lo==3||(lo==5)
            ylabel(ti{t},'FontSize',12);
        end
        hold off
        grid on;
        set(gca,'Layer','top','FontSize',12);
        %axis([0 ll-1 min(min(TE))-0.1 max(max(TE))+0.1])
        axis([0 ll-1 min(min(prctile(irf_chol(t,tt,:,:),5,4)))-.1 max(max(prctile(irf_chol(t,tt,:,:),95,4)))+.1]);
    end
end;
saveas(gcf,['NoisyBC_rev2' num2str(ff)],'png')

if length(varseries)>2
    ff=ff+1;
    figure(ff)
    lo=0;
    for t=[indGDP indC indI]
        for j=1:2,
            temp = squeeze(irf_chol(t,tt(j),1:ll,:));
            TEMP = irf_plot(t,tt(j),1:ll);
            TEMP2 = irf_plot_rev(t,tt(j),1:ll);
            lo=lo+1;
            subplot(3,2,lo);
            te2 = prctile(temp',[5])';
            te3 = prctile(temp',[95])';
            te4 = prctile(temp',[16])';
            te5 = prctile(temp',[84])';
            te_med = prctile(temp',[50])';
            b=fill([0:ll-1 ll-1:-1:0],[te2' flipud(te3)'],col1);set(b,'EdgeAlpha',0);hold on
            b=fill([0:ll-1 ll-1:-1:0],[te4' flipud(te5)'],col2);set(b,'EdgeAlpha',0);
            a = plot(0:ll-1,[squeeze(TEMP) squeeze(TEMP2) 0*te_med ]);
            set(a(1),'Color','k','LineWidth',2);
            set(a(2),'Color','r','LineWidth',2,'LineStyle','--');
            set(a(end),'Color','k');
            TE = [te2 te3 te4 te5];
            %c = plot(0:ll-1,TE,'k');
            if lo==1
                c=title('Surprise');set(c,'FontSize',12);
                ylabel(ti{t},'FontSize',12);
            elseif lo==2
                c=title('Signal');set(c,'FontSize',12);
            elseif (lo==3)||(lo==5)
                ylabel(ti{t},'FontSize',12);
            end
            hold off
            grid on;
            set(gca,'Layer','top','FontSize',12);
            %axis([0 ll-1 min(min(prctile(irfb(varseries(t),[ind1 ind],:,:),5,4)))-0.1 max(max(prctile(irfb(varseries(t),[ind1 ind],:,:),95,4)))+0.1 ])
            axis([0 ll-1 min(min(prctile(irf_chol(t,tt,:,:),5,4)))-.1 max(max(prctile(irf_chol(t,tt,:,:),95,4)))+.1]);
        end
    end;
    saveas(gcf,['NoisyBC_rev2' num2str(ff)],'png')
    
end

%% Long-run and noise

%irfsb = CI(varseries,:,:,:);
%se usi Kilian: irfsb = CI;

%irf_str = irfs(varseries,:,:);


ff=ff+1;

lo=0;
for t=[ind1 indControl ind]
    for j=1:2;
        temp = squeeze(irfsb(t,j,1:ll,:));
        
        TEMP = irf_str(t,j,1:ll);
        TEMP2 = irf_str_rev(t,j,1:ll);
        lo=lo+1;
        figure(ff)
        subplot(3,2,lo);
        
        %plot([squeeze(Cp(v_sel(t),1,1:40)) prctile(temp',[16 50 84])' prctile(temp',[5 95])']);
        
        te2 = prctile(temp',[5])';
        te3 = prctile(temp',[95])';
        te4 = prctile(temp',[16])';
        te5 = prctile(temp',[84])';
        te_med = prctile(temp',[50])';
        b=fill([0:ll-1 ll-1:-1:0],[te2' flipud(te3)'],col1);set(b,'EdgeAlpha',0);hold on
        b=fill([0:ll-1 ll-1:-1:0],[te4' flipud(te5)'],col2);set(b,'EdgeAlpha',0);
        a = plot(0:ll-1,[squeeze(TEMP) squeeze(TEMP2) 0*te_med ]);
            set(a(1),'Color','k','LineWidth',2);
            set(a(2),'Color','r','LineWidth',2,'LineStyle','--');
            set(a(end),'Color','k');
        TE = [te2 te3 te4 te5];
        %c = plot(0:ll-1,TE,'k');
        %a = plot([te_med]);set(a,'Color','k','LineWidth',2,'');
        if lo==1
            %c=title('News');set(c,'FontSize',12);
            c=title('Long-run');set(c,'FontSize',12);
            ylabel(ti{t},'FontSize',12);
        elseif lo==2
            c=title('Noise');set(c,'FontSize',12);
        elseif lo==3||lo==5
            ylabel(ti{t},'FontSize',12);
        end
        hold off
        grid on;
        set(gca,'Layer','top','FontSize',12);
        %axis(ax(t,:));
        %axis([0 ll-1 min(min(TE))-0.1 max(max(TE))+0.1 ])
        axis([0 ll-1 min(min(prctile(irfsb(t,1:2,:,:),5,4)))-.1 max(max(prctile(irfsb(t,1:2,:,:),95,4)))+.1]);
        %     xlabel('Quarters after shock','FontSize',10);
        %     ylabel('% Deviations','FontSize',10);
    end
end;
saveas(gcf,['NoisyBC_rev' num2str(ff)],'png')

if length(varseries)>2
    ff=ff+1;
    
    lo=0;
    for t=[indGDP indC indI]
        for j=1:2,
            temp = squeeze(irfsb(t,j,1:ll,:));
            
            TEMP = irf_str(t,j,1:ll);
            TEMP2 = irf_str_rev(t,j,1:ll);
            lo=lo+1;
            figure(ff);
            subplot(3,2,lo);
            
            %plot([squeeze(Cp(v_sel(t),1,1:40)) prctile(temp',[16 50 84])' prctile(temp',[5 95])']);
            
            te2 = prctile(temp',[5])';
            te3 = prctile(temp',[95])';
            te4 = prctile(temp',[16])';
            te5 = prctile(temp',[84])';
            te_med = prctile(temp',[50])';
            b=fill([0:ll-1 ll-1:-1:0],[te2' flipud(te3)'],col1);set(b,'EdgeAlpha',0);hold on
            b=fill([0:ll-1 ll-1:-1:0],[te4' flipud(te5)'],col2);set(b,'EdgeAlpha',0);
            a = plot(0:ll-1,[squeeze(TEMP) squeeze(TEMP2) 0*te_med ]);
            set(a(1),'Color','k','LineWidth',2);
            set(a(2),'Color','r','LineWidth',2,'LineStyle','--');
            set(a(end),'Color','k');
            TE = [te2 te3 te4 te5];
            %c = plot(0:ll-1,TE,'k');
            %a = plot([te_med]);set(a,'Color','k','LineWidth',2,'');
            if lo==1
                c=title('Long-run');set(c,'FontSize',12);
                %c=title('Tech');set(c,'FontSize',12);
                ylabel(ti{t},'FontSize',12);
            elseif lo==2
                c=title('Noise');set(c,'FontSize',12);
            elseif lo==3
                ylabel(ti{t},'FontSize',12);
            elseif lo==5
                ylabel(ti{t},'FontSize',12);
            end
            hold off
            grid on;
            set(gca,'Layer','top','FontSize',12);
            %axis(ax(t,:));
            %axis([0 ll-1 min(min(TE))-0.1 max(max(TE))+0.1 ])
            axis([0 ll-1 min(min(prctile(irfsb(t,1:2,:,:),5,4)))-.1 max(max(prctile(irfsb(t,1:2,:,:),95,4)))+.1]);
            %     xlabel('Quarters after shock','FontSize',10);
            %     ylabel('% Deviations','FontSize',10);
        end
    end;
    saveas(gcf,['NoisyBC_rev' num2str(ff)],'png')
    
end




