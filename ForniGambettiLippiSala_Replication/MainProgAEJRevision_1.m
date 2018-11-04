clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOAD DATA
rng('default')
plot_flag = 1;

[Data, codeData,tc] = ReadDataNews_EJ_data;

fa = 100;
Data(:,7) = fa*log(exp(Data(:,7)/fa)+exp(Data(:,13)/fa)); % inv
Data(:,11) = fa*log(exp(Data(:,12)/fa)+exp(Data(:,14)/fa));% cons

% prices in growth rates
pd = [31 33 34 35 36 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
temp2 = diff(Data(:,pd));
Data = Data(2:end,:);
Data(:,pd) = temp2;
Data = Data(4:end,:);

%fi = 8;
%out = filter(ones(1,fi)/fi,1,Data(:,93));
%Data(:,93) = out;disp('smoothed TFP')

Data(:,[93 94 97 98 99 100]) = cumsum(Data(:,[93 94 97 98 99 100]))/4;

pop = xlsread('population.xls');

if size(pop,2) ==2, 
    pop = pop(:,2);
end
pop_d = [1:18 20:21 26 95 96];
Data(:,pop_d) = fa*log(exp(Data(:,pop_d)/fa)./kron(pop(1:size(Data,1),1),ones(1,length(pop_d))));

Data(:,108) = 100*log(Data(:,108)./pop(1:size(Data,1))*1000000);

%Data = Data(fi:end,:);
Data = Data(3:end,:);

x = standardize(center(Data(:,1:end-3)));
pc = principalcomponents(x,15);


%% Parameters
kilian_flag = 1;
ll = 40;
npc = 0;
nrepli = 500;
lags = 4;
col2 = [.65 .65 .65];
col1 = [.85 .85 .85];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bivariate model
no = 1;

if no == 1;
    varseries = [108 83];
    codeData=ones(1,108);codeData(varseries(1)) = 0;
    %varseries = [93 103 ];
    
    ind1 = 1;
    ind = 2;
    [aic,bic,hqc] = aicbic(Data(:,varseries),4);
    [irfs,irfsb,irf,irfb,sigma_a,sigma_aboot,bidielle,sh,ssh,wr] = FAVARNewsNoiseChol(Data,...
        [Data(:,varseries) pc(:,1:npc)],codeData,lags,ll,nrepli,[varseries 10:10+npc-1],ind1,ind);
    % irfs: structural
    % irfsb: bootstrap replication
    % irf: reduced form (Cholesky)
    % irfb: bootstrap replication
    
    for k=[2 4]
        for j=[2 4 6]
            orto11(j/2,k/2) = ortotest(sh(1:end,1),diff(pc(lags:end,1:j)),k);% learning
            orto12(j/2,k/2) = ortotest(sh(1:end,2),diff(pc(lags:end,1:j)),k);% signal
            orto13(j/2,k/2) = ortotest(ssh(1:end-8,1),diff(pc(lags:end-8,1:j)),k);% long_run
            orto14(j/2,k/2) = ortotest(ssh(1:end-8,2),diff(pc(lags:end-8,1:j)),k);% noise
        end
    end
    
    ff = 0;
    
    %DoFiguresBC
    % close all
end
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Potential output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 93: TFP ut adjusted
% 94: TFP
%109: spread
%  1: GDP
% 11: C
%  7: I
%103: BC expected next 12 mths
%104: BC expcted next 5 yrs
%108: potential GDP
% 83: S&P
% 56: FFR

varseries = [108 56 83 1 11 7];

codeData = ones(1,109);codeData(varseries) = 0;

indGDP = find(varseries==1);
indC   = find(varseries==11);
indI   = find(varseries==7);
ind1   = find(varseries==108);% variabile colpita da news
indControl = find(varseries==56);% control
ind    = find(varseries==83);% signal

ti{ind1}='Potential GDP';ti{ind}='Stock Prices';ti{indGDP}='GDP';ti{indC}='Consumption';ti{indI}='Investment';ti{indControl}='FFR';

%[aic,bic,hqc] = aicbic(Data(:,varseries),4);
%lags=aic;
lags = 4;

nr = 500;
[irfs,irfsb,irf,irfb,sigma_a,sigma_aboot,bidielle,sh,ssh,wr] = FAVARNewsNoiseChol(Data,[Data(:,varseries) pc(:,1:npc)],codeData,lags,ll,nr,[varseries 15:15+npc-1],ind1,ind);
% irfs: structural
% irfsb: bootstrap replicationtab
% irf: reduced form (Cholesky)
% irfb: bootstrap replication
irfsbold=irfsb;
for k = 1:4
    for j = 1:10
        orto1(j,k) = ortotest(sh(1:end,1),diff(pc(lags:end,1:j)),k);% learning
        orto2(j,k) = ortotest(sh(1:end,2),diff(pc(lags:end,1:j)),k);% signal
        orto3(j,k) = ortotest(ssh(1:end-8,1),diff(pc(lags:end-8,1:j)),k);% long-run
        orto4(j,k) = ortotest(ssh(1:end-8,2),diff(pc(lags:end-8,1:j)),k);% noise
    end
end
%orto1([2 4 6], [2 4])'
%return
irfold = irf;
totvarnonstruc = squeeze(sum(cumsum(irfold.^2,3),2));

varu = squeeze(cumsum(irf(:,ind1,:).^2,3));
varalfa = squeeze(cumsum(irf(:,ind,:).^2,3));
varualfa = varu + varalfa;

tab_surprise = varu(varseries,[1 5 9 17  ll+1])./totvarnonstruc(varseries,[1 5 9 17 ll+1])
tab_signal = varalfa(varseries,[1 5 9 17  ll+1])./totvarnonstruc(varseries,[1 5 9 17 ll+1])

varnews0 = squeeze(cumsum(irfs(:,1,:).^2,3));
varnoise0 = squeeze(cumsum(irfs(:,2,:).^2,3));

varnewsnoise = varnews0 + varnoise0;
totvar = totvarnonstruc - varualfa + varnewsnoise;
varnews1 = varnews0./totvar;
varnoise1 = varnoise0./totvar;

tab_news = varnews1(varseries,[1 5 9 17  ll+1])
tab_noise = varnoise1(varseries,[1 5 9 17  ll+1])

if kilian_flag == 1;
    p = lags;
    X = Data(:,varseries);
    [A,SIGMA,U,V] = olsvarc(X,p);
    h = 40;
    [CI,CI_chol] = boot_luca(A,U,X,V,p,h,nr,codeData(varseries),ind,ind1);
    irf_chol = CI_chol;%(varseries,:,:,:);%se usi Kilian
    irfsb = CI;% %se usi Kilian: 
else
    CI_chol = irfb;
    CI = irfsb;
    irf_chol = CI_chol(varseries,:,:,:);
    irfsb = CI(varseries,:,:,:);
end

irf_plot = irf(varseries,:,:);

% structural form
irf_str = irfs(varseries,:,:);
irf_pot_bench = irf_str;% str IRFs

irfb_pot_bench = irfsb;% str bootstrap


save('IRF_pot.mat','irf_chol','irfsb','irf_plot','irf_str','irfb_pot_bench')

ff = 0;
HD_flag = 1;
DoFiguresBC_luca2_rev

%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TFP

varseries = [93 56 83 1 11 7];

codeData = ones(1,109);codeData(varseries) = 0;

indGDP = find(varseries==1);
indC   = find(varseries==11);
indI   = find(varseries==7);
ind1   = find(varseries==93);% variabile colpita da news
indControl = find(varseries==56);% control
ind    = find(varseries==83);% signal

ti{ind1}='TFP';ti{ind}='Stock Prices';ti{indGDP}='GDP';ti{indC}='Consumption';ti{indI}='Investment';ti{indControl}='FFR';

% [aic,bic,hqc] = aicbic(Data(:,varseries),4);
% lags=aic;
%lags = 4;

nr = 500;
[irfs,irfsb,irf,irfb,sigma_a,sigma_aboot,bidielle,sh,ssh,wr] = FAVARNewsNoiseChol(Data,[Data(:,varseries) pc(:,1:npc)],codeData,lags,ll,nr,[varseries 15:15+npc-1],ind1,ind);
% irfs: structural
% irfsb: bootstrap replicationtab
% irf: reduced form (Cholesky)
% irfb: bootstrap replication

for k = 1:4
    for j = 1:10
        orto5(j,k) = ortotest(sh(1:end,1),diff(pc(lags:end,1:j)),k);% learning
        orto6(j,k) = ortotest(sh(1:end,2),diff(pc(lags:end,1:j)),k);% signal
        orto7(j,k) = ortotest(ssh(1:end-8,1),diff(pc(lags:end-8,1:j)),k);% long_run
        orto8(j,k) = ortotest(ssh(1:end-8,2),diff(pc(lags:end-8,1:j)),k);% noise
    end
end

irfold = irf;
totvarnonstruc = squeeze(sum(cumsum(irfold.^2,3),2));

varu = squeeze(cumsum(irf(:,ind1,:).^2,3));
varalfa = squeeze(cumsum(irf(:,ind,:).^2,3));
varualfa = varu + varalfa;

tab_surprise = varu(varseries,[1 5 9 17  ll+1])./totvarnonstruc(varseries,[1 5 9 17  ll+1])
tab_signal = varalfa(varseries,[1 5 9 17  ll+1])./totvarnonstruc(varseries,[1 5 9 17  ll+1])

varnews0 = squeeze(cumsum(irfs(:,1,:).^2,3));
varnoise0 = squeeze(cumsum(irfs(:,2,:).^2,3));

varnewsnoise = varnews0 + varnoise0;
totvar = totvarnonstruc - varualfa + varnewsnoise;
varnews1 = varnews0./totvar;
varnoise1 = varnoise0./totvar;

tab_news = varnews1(varseries,[1 5 9 17  ll+1])
tab_noise = varnoise1(varseries,[1 5 9 17  ll+1])

if kilian_flag == 1;
    p = lags;
    X = Data(:,varseries);
    [A,SIGMA,U,V] = olsvarc(X,p);
    h = 40;
    [CI,CI_chol] = boot_luca(A,U,X,V,p,h,nr,codeData(varseries),ind,ind1);
    irf_chol = CI_chol;%(varseries,:,:,:);%se usi Kilian
    irfsb = CI;% %se usi Kilian: 
else
    CI_chol = irfb;
    CI = irfsb;
    irf_chol = CI_chol(varseries,:,:,:);
    irfsb = CI(varseries,:,:,:);
end

irf_plot = irf(varseries,:,:);

% structural form
irf_str = irfs(varseries,:,:);
irf_tfp_bench = irf_str;% str IRFs
irfb_tfp_bench = irfsb;% str bootstrap

ff = 9;
HD_flag = 2;
DoFiguresBC_luca2_rev

%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROBUSTNESS
%data = Data(51:end-1,:);
data = Data;

% potential output
varseries = [108 109 83 1 11 7];% per tornare indietro , metti 109 invece di 111
% a different  control: spread BAA-10YBond

codeData = ones(1,111);codeData(varseries) = 0;

indGDP = find(varseries==1);
indC   = find(varseries==11);
indI   = find(varseries==7);
ind1   = find(varseries==108);% variabile colpita da news
indControl = find(varseries==109);% control
ind    = find(varseries==83);% signal

ti{ind1}='Pot. GDP';ti{ind}='Stock Prices';ti{indGDP}='GDP';ti{indC}='Consumption';ti{indI}='Investment';ti{indControl}='FFR';

nr = 5;
[irfsR1,irfsbR1,irfR1,irfbR1,sigma_a,sigma_aboot,bidielle,sh,ssh,wr] = FAVARNewsNoiseChol(data,[data(:,varseries) pc(:,1:npc)],codeData,lags,ll,nr,[varseries 15:15+npc-1],ind1,ind);
irfold = irfR1;
totvarnonstruc = squeeze(sum(cumsum(irfold.^2,3),2));

varu = squeeze(cumsum(irfR1(:,ind1,:).^2,3));
varalfa = squeeze(cumsum(irfR1(:,ind,:).^2,3));
varualfa = varu + varalfa;

tab_surprise = varu(varseries,[1 5 9 17  ll+1])./totvarnonstruc(varseries,[1 5 9 17  ll+1])
tab_signal = varalfa(varseries,[1 5 9 17  ll+1])./totvarnonstruc(varseries,[1 5 9 17  ll+1])

varnews0 = squeeze(cumsum(irfsR1(:,1,:).^2,3));
varnoise0 = squeeze(cumsum(irfsR1(:,2,:).^2,3));

varnewsnoise = varnews0 + varnoise0;
totvar = totvarnonstruc - varualfa + varnewsnoise;
varnews1 = varnews0./totvar;
varnoise1 = varnoise0./totvar;

tab_news = varnews1(varseries,[1 5 9 17  ll+1])
tab_noise = varnoise1(varseries,[1 5 9 17  ll+1])
%[luR1,laR1,bcnoiseR1,bctotR1] = robustness_hist_dec(Data,lags,irfsR1,ssh,varseries,indGDP);
irfsR1 = irfsR1(varseries,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all
% potential output
varseries = [108 56  104 1 11 7];
% signal: 104 - BC expcted next 5 yrs

codeData = ones(1,111);codeData(varseries) = 0;

indGDP = find(varseries==1);
indC   = find(varseries==11);
indI   = find(varseries==7);
ind1   = find(varseries==108);% variabile colpita da news
indControl = find(varseries==56);% control
ind    = find(varseries==104);% signal

ti{ind1}='Pot. GDP';ti{ind}='Stock Prices';ti{indGDP}='GDP';ti{indC}='Consumption';ti{indI}='Investment';ti{indControl}='FFR';

nr = 5;
[irfsR2,irfsbR2,irfR2,irfbR2,sigma_a,sigma_aboot,bidielle,sh,ssh,wr] = FAVARNewsNoiseChol(Data,[Data(:,varseries) pc(:,1:npc)],codeData,lags,ll,nr,[varseries 15:15+npc-1],ind1,ind);
irfold = irfR2;
totvarnonstruc = squeeze(sum(cumsum(irfold.^2,3),2));

varu = squeeze(cumsum(irfR2(:,ind1,:).^2,3));
varalfa = squeeze(cumsum(irfR2(:,ind,:).^2,3));
varualfa = varu + varalfa;

tab_surprise = varu(varseries,[1 5 9 17  ll+1])./totvarnonstruc(varseries,[1 5 9 17  ll+1])
tab_signal = varalfa(varseries,[1 5 9 17  ll+1])./totvarnonstruc(varseries,[1 5 9 17  ll+1])

varnews0 = squeeze(cumsum(irfsR2(:,1,:).^2,3));
varnoise0 = squeeze(cumsum(irfsR2(:,2,:).^2,3));

varnewsnoise = varnews0 + varnoise0;
totvar = totvarnonstruc - varualfa + varnewsnoise;
varnews1 = varnews0./totvar;
varnoise1 = varnoise0./totvar;

tab_news = varnews1(varseries,[1 5 9 17  ll+1])
tab_noise = varnoise1(varseries,[1 5 9 17  ll+1])

%[luR2,laR2,bcnoiseR2,bctotR2] = robustness_hist_dec(Data,lags,irfsR2,ssh,varseries,indGDP);
irfsR2 = irfsR2([108 56 104 1 11 7],:,:);
%ff=100;
Robustness_luca2
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TFP + No news identification

varseries = [93 83 56 1 11 7];

codeData = ones(1,111);codeData(varseries) = 0;

indGDP = find(varseries==1);
indC   = find(varseries==11);
indI   = find(varseries==7);
ind1   = find(varseries==93);%93 variabile colpita da news
indControl = find(varseries==56);% control
ind    = find(varseries==83);% signal

ti{ind1}='TFP';ti{ind}='Stock Prices';ti{indGDP}='GDP';ti{indC}='Consumption';ti{indI}='Investment';ti{indControl}='FFR';

nr = 5;
[irfsR3,irfsbR3,irfR3,irfbR3,sigma_a,sigma_aboot,bidielle,sh,ssh,wr,irfchol] = FAVARNewsNoiseMax(Data,[Data(:,varseries) pc(:,1:npc)],codeData,lags,ll,nr,[varseries 15:15+npc-1],ind1,ind);
% irfs: structural
% irfsb: bootstrap replicationtab
% irf: reduced form (Cholesky)
% irfb: bootstrap replication

irfold = irfchol;
totvarnonstruc = squeeze(sum(cumsum(irfold.^2,3),2));

varu = squeeze(cumsum(irfR3(:,ind1,:).^2,3));
varalfa = squeeze(cumsum(irfR3(:,ind,:).^2,3));
varualfa = varu + varalfa;

tab_surprise = varu(varseries,[1 5 9 17  ll+1])./totvarnonstruc(varseries,[1 5 9 17  ll+1])
tab_signal = varalfa(varseries,[1 5 9 17  ll+1])./totvarnonstruc(varseries,[1 5 9 17  ll+1])

varnews0 = squeeze(cumsum(irfsR3(:,1,:).^2,3));
varnoise0 = squeeze(cumsum(irfsR3(:,2,:).^2,3));

varnewsnoise = varnews0 + varnoise0;
totvar = totvarnonstruc - varualfa + varnewsnoise;
varnews1 = varnews0./totvar;
varnoise1 = varnoise0./totvar;

tab_news = varnews1(varseries,[1 5 9 17  ll+1])
tab_noise = varnoise1(varseries,[1 5 9 17  ll+1])

irfsR3 = irfsR3(varseries([1 3 2 4 5 6]),:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varseries = [93 109 83  1 11 7];

codeData = ones(1,111);codeData(varseries) = 0;

indGDP = find(varseries==1);
indC   = find(varseries==11);
indI   = find(varseries==7);
ind1   = find(varseries==93);%93 variabile colpita da news
indControl = find(varseries==109);% control
ind    = find(varseries==83);% signal

ti{ind1}='TFP';ti{ind}='Stock Prices';ti{indGDP}='GDP';ti{indC}='Consumption';ti{indI}='Investment';ti{indControl}='FFR';

nr = 5;
[irfsR4,irfsbR4,irfR4,irfbR4,sigma_a,sigma_aboot,bidielle,sh,ssh,wr] = FAVARNewsNoiseChol(data,[data(:,varseries) pc(:,1:npc)],codeData,lags,ll,nr,[varseries 15:15+npc-1],ind1,ind);
% irfs: structural
% irfsb: bootstrap replicationtab
% irf: reduced form (Cholesky)
% irfb: bootstrap replication


irfold = irfR4;
totvarnonstruc = squeeze(sum(cumsum(irfold.^2,3),2));

varu = squeeze(cumsum(irfR4(:,ind1,:).^2,3));
varalfa = squeeze(cumsum(irfR4(:,ind,:).^2,3));
varualfa = varu + varalfa;

tab_surprise = varu(varseries,[1 5 9 17  ll+1])./totvarnonstruc(varseries,[1 5 9 17  ll+1])
tab_signal = varalfa(varseries,[1 5 9 17  ll+1])./totvarnonstruc(varseries,[1 5 9 17  ll+1])

varnews0 = squeeze(cumsum(irfsR4(:,1,:).^2,3));
varnoise0 = squeeze(cumsum(irfsR4(:,2,:).^2,3));

varnewsnoise = varnews0 + varnoise0;
totvar = totvarnonstruc - varualfa + varnewsnoise;
varnews1 = varnews0./totvar;
varnoise1 = varnoise0./totvar;

tab_news = varnews1(varseries,[1 5 9 17  ll+1])
tab_noise = varnoise1(varseries,[1 5 9 17  ll+1])


irfsR4 = irfsR4(varseries,:,:);


Robustness_mario3
