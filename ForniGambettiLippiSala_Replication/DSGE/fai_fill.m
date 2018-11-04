te2 = prctile(temp',[5])';
te3 = prctile(temp',[95])';
te4 = 0*prctile(temp',[16])';
te5 = 0*prctile(temp',[84])';
te_med = prctile(temp',[50])';
%te_med = mean(temp,2);

b = fill([0:ll-1 ll-1:-1:0],[te2' flipud(te3)'],col1);set(b,'EdgeAlpha',0);hold on
a = plot(0:39,[te_med]);set(a,'LineWidth',2,'LineStyle','--');axis([0 39 -Inf Inf]);hold on


