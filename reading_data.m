%*************************************************************************%
% Main                                                                    %
%                                                                         %
% NOTE let's avoid using 'end' to refer to a variable in dataset          %
%                                                                         %
% last change 6/6/2018                                                    %
%*************************************************************************%

clear
close all

%Data Info
% x_t|t     = first column of variables
% x_t+4|t   = fifth column of variables
% x_t|t-1   = second column of variables previous period
% x_t+4|t-1 = sixth column of variables previous period

filename = 'main_file';
sheet = 'Sheet1';
range = 'B1:AB300';
[dataset, var_names] = read_data2(filename, sheet, range);

for i = 1:size(dataset,2)
      eval([var_names{i} '= dataset(:,i);']);
end

%Building Zt
%Step 1 - Getting the forecasted growth rates
Delta_RGDP_t   = log(RGDP5_SPF) - log(RGDP1_SPF);
Delta_RDGP_t1  = log(RGDP6_SPF) - log(RGDP2_SPF);
%Delta_INDPROD_t = log(dataset(:,22)) - log(dataset(:,20));
%Delta_INDPROD_t1 = log(dataset(:,23)) - log(dataset(:,21));
%Investment is the sum between residential and non residential investment
%Delta_RINV_t = log(dataset(:,14) + dataset(:,18)) ...
%     - log(dataset(:,12) + dataset(:,16));
%Delta_RINV_t1 = log(dataset(:,15) + dataset(:,19)) ...
%      - log(dataset(:,13) + dataset(:,17));
%Step 2 - Revision in forecast growth rates
Z1 = Delta_RGDP_t(2:end) - Delta_RDGP_t1(1:end-1);
%Z2 = Delta_INDPROD_t(2:end) - Delta_INDPROD_t1(1:end-1);
%Z3 = Delta_RINV_t(2:end) - Delta_RINV_t1(1:end-1);

%Runniong OLS to obtain Ztilde
T              = size(Z1,1);
lag            = 4;
start          = lag;
lagged         = lag;
const          = ones(T-start-lagged,1);
X              = zeros(T-start-lagged,6+2*lag);
X(:,1)         = const;
X(:,2)         = MUNI1Y(start+1+1:end-lagged);
X(:,3)         = PDVMILY(start+1+1:end-lagged);
X(:,4)         = HAMILTON3YP(start+1+1:end-lagged);
X(:,5)         = RESID08(start+1+1:end-lagged);
X(:,6)         = TAXNARRATIVE(start+1+1:end-lagged);
for i = 1:2*lag+1
      X(:,6+i) = DTFP_UTIL(start-lag+i+1:end-lagged-lag+i-1);
end

Y                 = Z1(start+1:end-lagged);
[B, zhat, Ztilde] = quick_ols(Y,X);

Ztilde_graph = Ztilde + .05;
figure('Position',[100 100 1000 600])
figure(1)
area(Time(start+1+1:end-lagged),NBERDates(start+1+1:end-lagged),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
grid on
plot(Time(start+1+1:end-lagged),Ztilde_graph,'black-','Linewidth',3)
hold off
%xlim([12 252])
ylim([.03 .061])
%legend('NBER recessions','Weight on recession regime F(z)','Location','SouthOutside','Orientation','horizontal')
%legend('boxoff')


%*************************************************************************%
% 2nd stage - Smooth Transition Local Projections                         %
%                                                                         %
%*************************************************************************%
varlist = {'DTFP','Real GDP', 'Real Consumption', 'Unemployment'};
j = 3; %select variable

lags =8;
H = 20; %irfs horizon

%standardize Ztilde to get one std dev shock
Ztilde = Ztilde/std(Ztilde);
%stlp(y,x,u,fz(-1),lags,H); where y is the dep var, u is the shock, x are the controls


[IR_E, IR_R] = stlp(100*RealCons(start+1+1:end-lagged),0,Ztilde, ...
      ProbRecession(start:end-lagged-2),lags,H);


figure(2)
hold on
q = plot([1:H]',cumsum(IR_E), '-r', 'linewidth', 3);
h = plot([1:H]',cumsum(IR_R), '--b','linewidth', 3);
plot([1:H]', 0*[1:H]', ':k');
set(gca,'TickLabelInterpreter','latex')
title(varlist{j},'interpreter', 'latex', 'fontsize', 18);
xlabel('Quarter','interpreter','latex','fontsize',16);
ylabel('\% deviation from s.s.','interpreter','latex','fontsize',16);
l=legend([q h],{'Expansion','Recession'},'Location', 'SouthOutside','interpreter','latex');
set(l, 'box','on', 'FontSize',16,'Orientation','horizontal');
















