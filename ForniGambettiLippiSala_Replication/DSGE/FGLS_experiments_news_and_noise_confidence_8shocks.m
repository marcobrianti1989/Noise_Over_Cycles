clear Data

param = [0.4600    0.0400    0.7700    0.9500    0.9900    1.1900    1.4700    0.1800    0.3100    0.2900   11.0000    0.3500    0.5300    0.1900    2.0800    4.3300    3.4900    0.8800    0.8700 ...
         0         0    0.3000    0.4900    0.0500    0.9600    0.5500    1.5000    0.5000    1.2800    0.0250    0.9900];
param(27:28) = [1.013 .005];

simulation = 1;

if     simulation == 1 ; restrict_v = 0;SS = 7;param(32) = 0;               %noise shock has a positive variance + no additional shock;
elseif simulation == 2 ; restrict_v = 1;SS = 6;param(32) = 0;               %noise shock has a zero variance + no additional shock;
elseif simulation == 3 ; restrict_v = 0;SS = 7;param(32) = sqrt(2)*param(6);%noise shock has a positive variance + additional shock;
end

%    1         2     3     4       5      6         7      8       9        10     11      12      13      14      15       16     17       18      19      20      21     22       23       24     25      26       27       28               
%  [rho_d   rho_q  rho_p  rho_w  rho_g  sig_eps   sig_v   sig_p   sig_w   sig_g   sig_d   sig_q    h       alp     phi     chi     xi      cal     cal_w   iot     iot_w   mu_p    nu_p    mu_w    nu_w    rho_r   gam_pi   gam_y ]

obs = 500;
%[CUM_SER,SER,Z_SER,IRF,ZIRF,LIRF] = FGLS_4lags_generate(param,obs);
%Data = CUM_SER(:,1:7);

%x = standardize(center(Data));
%pc = principalcomponents(x,5);
    
%% Parameters
ll = 40;
npc = 0;
nrepli = 0;
lags = 8;

ind1 = 1; 
ind = 2;

%codeData = [0 0 0 0 0 0 1];
%nr = 1;
%[irfs,irfsb,irf,irfb,sigma_a,sigma_aboot,bidielle,sh,ssh,wr] = FAVARNewsNoiseChol(Data,[Data pc(:,1:npc)],codeData,lags,ll,nr,[1:7],ind1,ind);
% irfs: structural
% irfsb: bootstrap replication
% irf: reduced form (Cholesky)
% irfb: bootstrap replication

%for t=1:6;plot([squeeze(irfs(1+t,1,1:40)) squeeze(LIRF(t,1,1:40))]);title(t);pause;end
for t = 1:100
    %disp(t)
    %SS = 7;[CUM_SER,SER,Z_SER,IRF,ZIRF,LIRF] = FGLS_4lags_generate(param,obs);Data = CUM_SER(:,1:SS);codeData = [0 0 0 0 0 0 1];% with noise
    
    [CUM_SER,SER,Z_SER,IRF,ZIRF,LIRF] = FGLS_4lags_generate_restricted_8shocks(param,obs,restrict_v);
    Data = CUM_SER(:,1:SS);
    codeData = [0 0 0 0 0 0];% no noise
    save  Data
    %pause
    x = standardize(center(Data));
    pc = principalcomponents(x,5);
        
    load Conf
    if restrict_v == 0;
        Data = [Data(:,1) Conf Data(:,[2:4 6:7])];
    else;
        Data = [Data(:,1) Conf Data(:,[2:4 6])];
    end
    ind1 = 1;ind = 2;
    nr = 1;
    [irfs,irfsb,irf,irfb,sigma_a,sigma_aboot,bidielle,sh,ssh,wr] = FAVARNewsNoiseChol(Data,[Data pc(:,1:npc)],codeData,lags,ll,nr,[1:SS],ind1,ind);
    
    IRFS(:,:,:,t) = irfs;
    
    %[m,C] = modulrootsNdim(im,series)
end
%close all
%return

if SS == 7;
    ser{1}='a';ser{2}='Conf';ser{3}='Y';ser{4}='C';ser{5}='I';ser{6}='Hours';ser{7}='Inflation';
else
    ser{1}='a';ser{2}='Conf';ser{3}='Y';ser{4}='C';ser{5}='I';ser{6}='Hours';
end

flag_fig_sim = 0;

if flag_fig_sim == 1;
    if SS == 7;
        h = figure(1)
        lo = 0;
        st = 0;
        co='Color';r='r';ls='LineStyle';lw='LineWidth';tr='--';wi=1.5;sa='set(a(2:3),co,r,ls,tr,lw,wi)';
        for sh = 1:2
            subplot(2,SS,1+st);a = plot([squeeze(LIRF(1,sh,1:40)) prctile(squeeze(IRFS(3,sh,1:40,:))',[2.5 97.5])']);set(a(1),'LineWidth',2);eval(sa);title(ser{3});grid;if sh==1;ylabel('Real');else;ylabel('Noise');end
            subplot(2,SS,2+st);a = plot([squeeze(LIRF(2,sh,1:40)) prctile(squeeze(IRFS(4,sh,1:40,:))',[2.5 97.5])']);set(a(1),'LineWidth',2);eval(sa);title(ser{4});grid;
            subplot(2,SS,3+st);a = plot([squeeze(LIRF(3,sh,1:40)) prctile(squeeze(IRFS(5,sh,1:40,:))',[2.5 97.5])']);set(a(1),'LineWidth',2);eval(sa);title(ser{5});grid;
            subplot(2,SS,4+st);a = plot([squeeze(LIRF(5,sh,1:40)) prctile(squeeze(IRFS(6,sh,1:40,:))',[2.5 97.5])']);set(a(1),'LineWidth',2);eval(sa);title(ser{6});grid;
            subplot(2,SS,5+st);a = plot([squeeze(LIRF(6,sh,1:40)) prctile(squeeze(IRFS(7,sh,1:40,:))',[2.5 97.5])']);set(a(1),'LineWidth',2);eval(sa);title(ser{7});grid;
            subplot(2,SS,6+st);a = plot([squeeze(IRF(20,sh,1:40)) prctile(squeeze(IRFS(1,sh,1:40,:))',[2.5 97.5])']);set(a(1),'LineWidth',2);eval(sa);title(ser{1});grid;% irf di a_t
            subplot(2,SS,7+st);a = plot([prctile(squeeze(IRFS(2,sh,1:40,:))',[2.5 97.5])']);set(a(1:2),'Color','r','LineStyle','--','LineWidth',wi);title(ser{2});grid;% irf di a_t
            st=SS;
        end
        %saveas(h,'model_noise.png');
    end
    
    if SS == 6
        h=figure(3);
        lo = 0;
        st = 0;
        for sh = 1:2
            subplot(2,SS,1+st);a = plot([squeeze(LIRF(1,sh,1:40)) prctile(squeeze(IRFS(3,sh,1:40,:))',[2.5 97.5])']);set(a(1),'LineWidth',2);set(a(2:3),'Color','r','LineStyle','--','LineWidth',1.5);title(ser{3});grid;if sh==1;ylabel('Real');else;ylabel('Noise');end
            subplot(2,SS,2+st);a = plot([squeeze(LIRF(2,sh,1:40)) prctile(squeeze(IRFS(4,sh,1:40,:))',[2.5 97.5])']);set(a(1),'LineWidth',2);set(a(2:3),'Color','r','LineStyle','--','LineWidth',1.5);title(ser{4});grid;
            subplot(2,SS,3+st);a = plot([squeeze(LIRF(3,sh,1:40)) prctile(squeeze(IRFS(5,sh,1:40,:))',[2.5 97.5])']);set(a(1),'LineWidth',2);set(a(2:3),'Color','r','LineStyle','--','LineWidth',1.5);title(ser{5});grid;
            subplot(2,SS,4+st);a = plot([squeeze(LIRF(5,sh,1:40)) prctile(squeeze(IRFS(6,sh,1:40,:))',[2.5 97.5])']);set(a(1),'LineWidth',2);set(a(2:3),'Color','r','LineStyle','--','LineWidth',1.5);title(ser{6});grid;
            subplot(2,SS,5+st);a = plot([squeeze(IRF(20,sh,1:40)) prctile(squeeze(IRFS(1,sh,1:40,:))',[2.5 97.5])']);set(a(1),'LineWidth',2);set(a(2:3),'Color','r','LineStyle','--','LineWidth',1.5);title(ser{1});grid;% irf di a_t
            %subplot(2,SS,6+st);a = plot([prctile(squeeze(IRFS(2,sh,1:40,:))',[2.5 97.5])']);set(a(1:2),'Color','r','LineStyle','--','LineWidth',1.5);title(ser{2});grid;% irf di a_t
            st=SS;
        end
        %saveas(h,'model_no_noise.png');
    end
    
end

%close all


if SS == 7;
    ser{1}='a';ser{2}='Conf';ser{3}='Y';ser{4}='C';ser{5}='I';ser{6}='Hours';ser{7}='Inflation';
else
    ser{1}='a';ser{2}='Conf';ser{3}='Y';ser{4}='C';ser{5}='I';ser{6}='Hours';
end

TO = 0;
col1 = [.85 .85 .85];
%figure
if SS == 7;
    h = figure(200)
    lo = 0;
    st = 0;
    for sh = 1:2
        %temp = prctile(squeeze(IRFS(3,sh,1:40,:))',[2.5 97.5])';
        temp = squeeze(IRFS(3,sh,1:40,:));
        subplot(SS-1-TO,2,1+st);fai_fill;a = plot(0:39,[squeeze(LIRF(1,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');grid;%axis([0 39 -1 3]);
        if sh==1;ylabel(ser{3});end;
        if sh==1;title('Long-run');else;title('Noise');end;        
        
        temp = squeeze(IRFS(4,sh,1:40,:));
        subplot(SS-1-TO,2,3+st);fai_fill;a = plot(0:39,[squeeze(LIRF(2,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');
        if sh==1;ylabel(ser{4});end;grid;%axis([0 39 -1 3]);
        
        temp = squeeze(IRFS(5,sh,1:40,:));
        subplot(SS-1-TO,2,5+st);fai_fill;a = plot(0:39,[squeeze(LIRF(3,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');
        if sh==1;ylabel(ser{5});end;grid;%axis([0 39  -2 6]);
        
        temp = squeeze(IRFS(6,sh,1:40,:));
        subplot(SS-1,2,7+st);fai_fill;a = plot(0:39,[squeeze(LIRF(5,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');
        if sh==1;ylabel(ser{6});end;grid;%axis([0 39  -1 1.5]);
        
        temp = squeeze(IRFS(7,sh,1:40,:));
        subplot(SS-1,2,9+st);fai_fill;a = plot(0:39,[squeeze(LIRF(6,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');
        if sh==1;ylabel(ser{7});end;grid;%axis([0 39  -.2 .1]);
                        
        temp = squeeze(IRFS(1,sh,1:40,:));
        subplot(SS-1,2,11+st);fai_fill;a = plot(0:39,[squeeze(IRF(20,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');
        %subplot(SS-1-TO,2,7+st);fai_fill;a = plot(0:39,[squeeze(IRF(20,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');
        if sh==1;ylabel(ser{1});end;grid;%axis([0 39  -1 2.5]);
        
        st = 1;
        
    end
    if simulation == 1;saveas(h,'model_noise_rev.png');elseif simulation == 3;saveas(h,'model_noise_8_shocks_rev.png');close;end
end


if SS == 6
    h=figure(400);
    lo = 0;
    st = 0;
    co='Color';r='r';ls='LineStyle';lw='LineWidth';tr='--';sa='set(a(2:3),co,r,ls,tr,lw,1.5)';
    for sh = 1:2
        temp = squeeze(IRFS(3,sh,1:40,:));
        subplot(SS-1,2,1+st);fai_fill;a = plot(0:39,[squeeze(LIRF(1,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');if sh==1;ylabel(ser{3});end;grid;%axis([0 39 -1 3]);
        if sh==1;title('Long-run');else;title('Noise');end
        
        temp = squeeze(IRFS(4,sh,1:40,:));
        subplot(SS-1,2,3+st);fai_fill;a = plot(0:39,[squeeze(LIRF(2,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');if sh==1;ylabel(ser{4});end;grid;%axis([0 39 -1 3]);
        
        temp = squeeze(IRFS(5,sh,1:40,:));
        subplot(SS-1,2,5+st);fai_fill;a = plot(0:39,[squeeze(LIRF(3,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');if sh==1;ylabel(ser{5});end;grid;%axis([0 39 -2 6]);
        
        temp = squeeze(IRFS(6,sh,1:40,:));
        subplot(SS-1,2,7+st);fai_fill;a = plot(0:39,[squeeze(LIRF(5,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');if sh==1;ylabel(ser{6});end;grid;%axis([0 39 -1 1.5]);
        
        %temp = prctile(squeeze(IRFS(7,sh,1:40,:))',[2.5 97.5])';
        %subplot(SS-1,2,9+st);fai_fill;a = plot(0:39,[squeeze(LIRF(6,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');if sh==1;ylabel(ser{3});end;grid;
                        
        temp = squeeze(IRFS(1,sh,1:40,:));
        subplot(SS-1,2,9+st);fai_fill;a = plot(0:39,[squeeze(IRF(20,sh,1:40))]);set(a(1),'LineWidth',2,'Color','k');if sh==1;ylabel(ser{1});end;grid;%axis([0 39 -1 2.5]);

        st = 1;
    end
    if simulation == 2;saveas(h,'model_no_noise_rev.png');close;end
end
%close all
