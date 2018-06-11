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
range = 'B1:AA290';
[dataset, var_names] = read_data2(filename, sheet, range);

%Building Zt 
%Step 1 - Getting the forecasted growth rates
Delta_RGDP_t   = log(dataset(:,10)) - log(dataset(:,8));
Delta_RDGP_t1  = log(dataset(:,11)) - log(dataset(:,9));
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
dataset   = dataset(2:end,:);

%Runniong OLS to obtain Ztilde
[T, nvar] = size(dataset);
lag = 11;
start     = lag+1;
lagged    = lag;
%lagged    = 0;
const     = ones(T-start-lagged+1,1);
% X         = [const dataset(start:end-lagged,2:7) ...
%     dataset(start-1:end-lagged-1,4)  ...
%     dataset(start-3:end-lagged-3,4) ...
%     dataset(start-4:end-lagged-4,4)  ...
%     dataset(start-5:end-lagged-5,4)  ...
%     dataset(start-6:end-lagged-6,4)  ...
%     dataset(start-7:end-lagged-7,4) ...
%     dataset(start-8:end-lagged-8,4) ...
%     ];

X         = [const dataset(start:end-lagged,2:7) ...
dataset(start+2:end-lagged+2,4) ...
      dataset(start+1:end-lagged+1,4) dataset(start-2:end-lagged-2,4) ...
      dataset(start-1:end-lagged-1,4) dataset(start+3:end-lagged+3,4) ...
      dataset(start-3:end-lagged-3,4) dataset(start+4:end-lagged+4,4) ...
      dataset(start-4:end-lagged-4,4) dataset(start+5:end-lagged+5,4) ...
      dataset(start-5:end-lagged-5,4) dataset(start+6:end-lagged+6,4) ...
     dataset(start-6:end-lagged-6,4) ...
dataset(start+7:end-lagged+7,4) ...
      dataset(start-7:end-lagged-7,4) dataset(start+8:end-lagged+8,4) ...
      dataset(start-8:end-lagged-8,4) dataset(start+9:end-lagged+9,4) ...
      dataset(start-9:end-lagged-9,4) dataset(start+10:end-lagged+10,4) ...
      dataset(start-10:end-lagged-10,4) dataset(start+11:end-lagged+11,4) ...
      dataset(start-11:end-lagged-11,4)
];
Y         = Z1(start:end-lagged);
[B, zhat, Ztilde] = quick_ols(Y,X);

% figure(1)
% hold on
% grid on
% plot(dataset(3:end-2,1),Ztilde)
% hold off

% corr(Ztilde(1:end),dataset(start:end-lagged,4))
for fo = 1:20
    cori(fo) = corr(Ztilde(1:end-fo),dataset(start+fo:end-lagged,4));
end
% plot(cori)

% figure
% plot(dataset(start:end-lagged,1),Ztilde)

Ztilde_graph = Ztilde + .05;
figure('Position',[100 100 1000 600])
figure(1) 
area(dataset(start:end-lagged,1),dataset(start:end-lagged,24),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
grid on
plot(dataset(start:end-lagged,1),Ztilde_graph,'black-','Linewidth',3)
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
j = 2; %select variable

lags =8; 
H = 20; %irfs horizon

%standardize Ztilde to get one std dev shock
Ztilde = Ztilde/std(Ztilde);
%stlp(y,x,u,fz(-1),lags,H); where y is the dep var, u is the shock, x are the controls


[IR_E, IR_R] = stlp(100*dataset(start:end-lagged,26),0,Ztilde, dataset(start-1:end-lagged-1,25),lags,H);

%100*dataset(start:end-lagged,[2 3 5 6 7])

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
















