clear
%close all
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


fi = 32;
out = filter(ones(1,fi)/fi,1,(Data(:,1)));% applied to 100*log output
me = mean([(Data(fi+1:end,1)) 100*log(Data(fi+1:end,108)) out(fi+1:end)]);
%OUT = (out(fi+1:end)-me(3)+me(2));
OUT = out(1:end);
%plot([(Data(fi+1:end,1)) 100*log(Data(fi+1:end,108)) out(fi+1:end)-me(3)+me(2)]);% both in logs
%plot([(Data(fi+1:end,1))-100*log(Data(fi+1:end,108)) (Data(fi+1:end,1))-(out(fi+1:end)-me(3)+me(2))]);% cyclical output


Data(:,[93 94 97 98 99 100]) = cumsum(Data(:,[93 94 97 98 99 100]))/4;

pop = xlsread('population.xls');

if size(pop,2) ==2, 
    pop = pop(:,2);
end
pop_d = [1:18 20:21 26 95 96];
Data(:,pop_d) = fa*log(exp(Data(:,pop_d)/fa)./kron(pop(1:size(Data,1),1),ones(1,length(pop_d))));


CBO_pot = 100*log(Data(:,108)./pop(1:size(Data,1)));% cbp_pot/population and logged
OUT2 = 100*log((exp(OUT/100))./pop(1:size(Data,1)));
me = mean([(Data(fi:end,1)) CBO_pot(fi:end) OUT2(fi:end)]);


figure(1000);plot(([CBO_pot(fi:end,1) OUT2(fi:end)-me(3)+me(2) Data(fi:end,1)]));grid
figure(1001);plot(([Data(fi:end,1)-CBO_pot(fi:end,1) Data(fi:end,1)-(OUT2(fi:end)-me(3)+me(2))]));grid
pause(2);close all
Data(:,108) = OUT2-me(3)+me(2);

Data = Data(fi:end,:);
%Data = Data(3:end,:);

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

if no == 0;
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
    irf_chol_rev = CI_chol;%(varseries,:,:,:);%se usi Kilian
    irfsb_rev = CI;% %se usi Kilian: 
else
    CI_chol = irfb;
    CI = irfsb;
    irf_chol_rev = CI_chol(varseries,:,:,:);% bootstrap choleski
    irfsb_rev = CI(varseries,:,:,:);% bootstrap structural
end

irf_plot_rev = irf(varseries,:,:);% point est choleski

% structural form
irf_str_rev = irfs(varseries,:,:);% 
irf_pot_bench__rev = irf_str_rev;% point est structural IRFs

irfb_pot_bench_rev = irfsb;% structural bootstrap

load('IRF_pot.mat','irf_chol','irfsb','irf_plot','irf_str','irfb_pot_bench');

ff = 0;
HD_flag = 1;
DoFiguresBC_luca2_rev2

