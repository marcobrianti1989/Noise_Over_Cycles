%irf_pot_bench - str. IRFs
%irfb_pot_bench - str. bootstrap

ff=ff+1;
figure(ff)
lo=0;
for t=[ind1 indControl ind]
    for j=1:2;
        temp = squeeze(irfb_pot_bench(t,j,1:ll,:));
        if t==2
            TEMP = [squeeze(irf_pot_bench(t,j,1:ll)) 0*squeeze(irfsR1(t,j,1:ll)) squeeze(irfsR2(t,j,1:ll))];
        elseif t==3
            TEMP = [squeeze(irf_pot_bench(t,j,1:ll)) squeeze(irfsR1(t,j,1:ll)) 0*squeeze(irfsR2(t,j,1:ll))];
        else
        TEMP = [squeeze(irf_pot_bench(t,j,1:ll)) squeeze(irfsR1(t,j,1:ll)) squeeze(irfsR2(t,j,1:ll))];
        end
        lo=lo+1;
        
        subplot(3,2,lo);
        
        %plot([squeeze(Cp(v_sel(t),1,1:40)) prctile(temp',[16 50 84])' prctile(temp',[5 95])']);
        
        te2 = prctile(temp',[5])';
        te3 = prctile(temp',[95])';
        te4 = 0*prctile(temp',[16])';
        te5 = 0*prctile(temp',[84])';
        te_med = prctile(temp',[50])';
        b=fill([0:ll-1 ll-1:-1:0],[te2' flipud(te3)'],col1);set(b,'EdgeAlpha',0);hold on
        b=fill([0:ll-1 ll-1:-1:0],[te4' flipud(te5)'],col2);set(b,'EdgeAlpha',0);
        a = plot(0:ll-1,[squeeze(TEMP) 0*te_med ]);
        set(a(1),'Color','k','LineWidth',2);%if lo==2;legend(a(1:3),'Benchmark specification I','Stock prices with spread','Business cond. 5 yrs');end
        if t==2
            set(a(2),'Color','k');
        else
        set(a(2),'LineWidth',2,'LineStyle','--');
        end
        if t==3
            set(a(3),'Color','k');
        else
        set(a(3),'LineWidth',1,'LineStyle','-','Marker','x');
        end
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
        axis([0 ll-1 min(min(prctile(irfb_pot_bench(t,1:2,:,:),5,4)))-.1 max(max(prctile(irfb_pot_bench(t,1:2,:,:),95,4)))+.1]);
        %     xlabel('Quarters after shock','FontSize',10);
        %     ylabel('% Deviations','FontSize',10);
    end
end;
saveas(gcf,['NoisyBC' num2str(ff)],'png')

ff = ff+1; 
figure(ff);
if length(varseries)>2
       
    lo = 0;
    for t = [indGDP indC indI]
        for j = 1:2,
            temp = squeeze(irfb_pot_bench(t,j,1:ll,:));
            
            TEMP = [squeeze(irf_pot_bench(t,j,1:ll)) squeeze(irfsR1(t,j,1:ll)) squeeze(irfsR2(t,j,1:ll))];            
            lo = lo+1;
            
            subplot(3,2,lo);
            
            %plot([squeeze(Cp(v_sel(t),1,1:40)) prctile(temp',[16 50 84])' prctile(temp',[5 95])']);
            
            te2 = prctile(temp',[5])';
            te3 = prctile(temp',[95])';
            te4 = 0*prctile(temp',[16])';
            te5 = 0*prctile(temp',[84])';
            te_med = prctile(temp',[50])';
            b=fill([0:ll-1 ll-1:-1:0],[te2' flipud(te3)'],col1);set(b,'EdgeAlpha',0);hold on
            b=fill([0:ll-1 ll-1:-1:0],[te4' flipud(te5)'],col2);set(b,'EdgeAlpha',0);
            a = plot(0:ll-1,[squeeze(TEMP) 0*te_med ]);
            set(a(1),'Color','k','LineWidth',2);%if lo==2;legend(a(1:3),'Benchmark specification I','Stock prices with spread','Business cond. 5 yrs');end
            set(a(2),'LineWidth',2,'LineStyle','--');
            set(a(3),'LineWidth',1,'LineStyle','-','Marker','x');
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
            axis([0 ll-1 min(min(prctile(irfb_pot_bench(t,1:2,:,:),5,4)))-.1 max(max(prctile(irfb_pot_bench(t,1:2,:,:),95,4)))+.1]);
            %     xlabel('Quarters after shock','FontSize',10);
            %     ylabel('% Deviations','FontSize',10);
        end
    end;
    saveas(gcf,['NoisyBC' num2str(ff)],'png')
    
end




return

tim = 1960.25+0.25*lags+0.25*4:.25:2010.75;
start = 128;
stop = 4;

timespan = 70;
bcnoise = cfilter(noisecomp,6,32);
bctot = cfilter(datacomp,6,32);

ff = ff+1;
figure(ff)
subplot(2,1,1),plot(tim(end-timespan:end-stop),lu(end-timespan:end-stop),'k', ...
    tim(end-timespan:end-stop),lu(end-timespan:end-stop)-la(end-timespan:end-stop),'--r','Linewidth',1.5 ),grid on

subplot(2,1,2),,plot(tim(end-timespan:end-stop),bcnoise(end-timespan:end-stop),'--r',...
    tim(end-timespan:end-stop), bctot(end-timespan:end-stop),'k','Linewidth',1.5 ),grid on
saveas(gcf,['NoisyBC' num2str(ff)],'png')
close

