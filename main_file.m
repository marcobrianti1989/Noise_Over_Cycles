%*************************************************************************%
% Main                                                                    %
%                                                                         %
% NOTE describe variables (especially SHOCKS) in dataset        %
%                                                                         %
% last change 6/12/2018                                                    %
%*************************************************************************%

clear
close all

%Data Info
% x_t|t     = first column of variables
% x_t+4|t   = fifth column of variables
% x_t|t-1   = second column of variables previous period
% x_t+4|t-1 = sixth column of variables previous period

%Read main dataset
filename = 'main_file';
sheet = 'Sheet1';
range = 'B1:CC300';
do_truncation = 0; %Do not truncate data. You will have many NaN
[dataset, var_names] = read_data2(filename, sheet, range, do_truncation);
dataset = real(dataset);
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

%Read dataset_PC for PC analysis
filename_PC = 'Dataset_test_PC';
sheet_PC = 'Quarterly';
range_PC = 'B2:DA300';
do_truncation_PC = 1; %Do truncate data.
[dataset_PC, var_names_PC] = read_data2(filename_PC, sheet_PC, range_PC, do_truncation_PC);
dataset_PC = real(dataset_PC);
date_start_PC = dataset_PC(1,1);
dataset_PC = dataset_PC(:,2:end); %Removing time before PC analysis
PC = get_principal_components(dataset_PC);
pc = nan(size(dataset,1),size(dataset_PC,2));
loc_time_PC = find(Time == date_start_PC);
pc(loc_time_PC:loc_time_PC+size(PC,1)-1,:) = PC;


%Building Zt
%Step 1 - Getting the forecasted growth rates
%Real GDP
Delta_RGDP_t        = log(RGDP5_SPF) - log(RGDP1_SPF);
Delta_RDGP_t1       = log(RGDP6_SPF) - log(RGDP2_SPF);
%Industrial Production
Delta_INDPROD_t     = log(dataset(:,22)) - log(dataset(:,20));
Delta_INDPROD_t1    = log(dataset(:,23)) - log(dataset(:,21));
%Investment is the sum between residential and non residential investment
Delta_RINV_t        = log(dataset(:,14) + dataset(:,18)) ...
      - log(dataset(:,12) + dataset(:,16));
Delta_RINV_t1       = log(dataset(:,15) + dataset(:,19)) ...
      - log(dataset(:,13) + dataset(:,17));
%Step 2 - Revision in forecast growth rates
Z1 = Delta_RGDP_t(2:end) - Delta_RDGP_t1(1:end-1);
Z2 = Delta_INDPROD_t(2:end) - Delta_INDPROD_t1(1:end-1);
Z3 = Delta_RINV_t(2:end) - Delta_RINV_t1(1:end-1);

Z  = [Z1 Z2 Z3];

for iii = 1:1
      
      if iii == 4
            disp('Industrial Production does not work.')
      else
            ZZ = Z(:,iii);
            
            %First Stage - Building Ztilde
            lag_tfp        = 8; %number of lags of TFP - cannot be zero since 1 include current TFP
            lead_tfp       = 16; %number of leads of TFP
            lag            = 2;  %number of lags of control variables (other structural shocks)
            mpc            = 2; %max number of principal components
            threshold      = -1/eps;
            loc_start      = find(ZZ > threshold, 1) + lag;
            loc_end        = find(isnan(MUNI1Y(loc_start+1:end)),1);
            loc_end        = loc_start + loc_end - 1-2;
            ZZ             = ZZ(loc_start:loc_end-1);
%             for it = 1:10
%                   cori(it) = corr(ZZ,DTFP_UTIL(loc_start+1+it:loc_end+it))
%             end
%             plot(cori)
           
            
            %Runniong OLS to obtain Ztilde
            T              = size(ZZ,1);
            const          = ones(T,1);
            %X             = zeros(T,6 + 2*lag_tfp + 1);
            X              = const;
            
            DATA = [MUNI1Y, PDVMILY, HAMILTON3YP, RESID08, TAXNARRATIVE];
            trend = 1:1:length(X);
            X = [X,trend' , DATA(loc_start+1:loc_end,:)];
            SX = size(X,2);
            for i = 1:lag_tfp %When i = 1 TFP is contemporaneous
                  %X(:,end+1) = DTFP_UTIL(loc_start+2-i:loc_end-i+1);
                  X(:,end+1) = DTFP_UTIL(loc_start+2-i:loc_end-i+1).*ProbRecession(loc_start:loc_end-1);
            end
            SX = size(X,2);
            for i = 1:lead_tfp
                  %X(:,end+1) = DTFP_UTIL(loc_start+1+i:loc_end+i);
                  X(:,end+1) = DTFP_UTIL(loc_start+1+i:loc_end+i).*ProbRecession(loc_start:loc_end-1);
            end
            for l = 1:lag
                  X = [X DATA(loc_start+1-l:loc_end-l,:) pc(loc_start+1-l:loc_end-l,1:mpc)];
            end
            Y = ZZ;
%             Y = 0.1*Mich1Y(loc_start+1:loc_end);
%             Y = 0.1*SP500(loc_start+1:loc_end);
            [B, zhat, Ztilde] = quick_ols(Y,X);
                        for it = 1:10
                  cori(it) = corr(Ztilde.*ProbRecession(loc_start:loc_end-1),DTFP_UTIL(loc_start+1+it:loc_end+it))
            end
            plot(cori)
            return
            
            Ztilde_graph = Ztilde + .05;
            figure('Position',[100 100 1000 600])
            figure(iii)
            area(Time(loc_start+1:loc_end),NBERDates(loc_start+1:loc_end),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
            hold on
            grid on
            plot(Time(loc_start+1:loc_end),Ztilde_graph,'black-','Linewidth',3)
            hold off
            %xlim([12 252])
            ylim([.03 .061])
            %legend('NBER recessions','Weight on recession regime F(z)','Location','SouthOutside','Orientation','horizontal')
            %legend('boxoff')
            
            
            %*************************************************************************%
            % 2nd stage - Smooth Transition Local Projections                         %
            %                                                                         %
            %*************************************************************************%
            
            % STLP
            varlist = {'TFP','Real GDP', 'Real Consumption', 'Unemployment Rate','Real Wage','Hours','CPI','Real Investment'};
            lags =1;
            H = 20; %irfs horizon
            mpc = 2;
            
            %standardize Ztilde to get one std dev shock
            Ztilde = Ztilde/std(Ztilde);
            %stlp(y,x,u,fz(-1),lags,H); where y is the dep var, u is the shock, x are the controls
            
            [IR_E_C, IR_R_C, IR_L_C] = stlp(100*RealCons(loc_start+1:loc_end-2),...
                  [pc(loc_start+1:loc_end-2,1:mpc)] ,Ztilde(1:end-2), ...
                  ProbRecession(loc_start:loc_end-1-2),lags,H,DTFP_UTIL(loc_start+1:loc_end-2));
            
            [IR_E_G, IR_R_G, IR_L_G] = stlp(100*RealGDP(loc_start+1:loc_end-2),[pc(loc_start+1:loc_end-2,1:mpc)],Ztilde(1:end-2), ...
                  ProbRecession(loc_start:loc_end-1-2),lags,H,DTFP_UTIL(loc_start+1:loc_end-2));
            
            [IR_E_T, IR_R_T, IR_L_T] = stlp(100*DTFP_UTIL(loc_start+1:loc_end-2),[pc(loc_start+1:loc_end-2,1:mpc)],Ztilde(1:end-2), ...
                  ProbRecession(loc_start:loc_end-1-2),lags,H,DTFP_UTIL(loc_start+1:loc_end-2));
            
            [IR_E_U, IR_R_U, IR_L_U] = stlp(UnempRate(loc_start+1:loc_end-2),[pc(loc_start+1:loc_end-2,1:mpc)],Ztilde(1:end-2), ...
                  ProbRecession(loc_start:loc_end-1-2),lags,H,DTFP_UTIL(loc_start+1:loc_end-2));
            
            [IR_E_W, IR_R_W, IR_L_W] = stlp(100*RealWage(loc_start+1:loc_end-2),[pc(loc_start+1:loc_end-2,1:mpc)],Ztilde(1:end-2), ...
                  ProbRecession(loc_start:loc_end-1-2),lags,H,DTFP_UTIL(loc_start+1:loc_end-2));
            
            [IR_E_H, IR_R_H, IR_L_H] = stlp(100*Hours(loc_start+1:loc_end-2),[pc(loc_start+1:loc_end-2,1:mpc)],Ztilde(1:end-2), ...
                  ProbRecession(loc_start:loc_end-1-2),lags,H,DTFP_UTIL(loc_start+1:loc_end-2));
            
            [IR_E_CPI, IR_R_CPI, IR_L_CPI] = stlp(100*CPI(loc_start+1:loc_end-2),[pc(loc_start+1:loc_end-2,1:mpc)],Ztilde(1:end-2), ...
                  ProbRecession(loc_start:loc_end-1-2),lags,H,DTFP_UTIL(loc_start+1:loc_end-2));
            
            [IR_E_INV, IR_R_INV, IR_L_INV] = stlp(100*RealInvestment(loc_start+1:loc_end-2),[pc(loc_start+1:loc_end-2,1:mpc)],Ztilde(1:end-2), ...
                  ProbRecession(loc_start:loc_end-1-2),lags,H,DTFP_UTIL(loc_start+1:loc_end-2));
            
            [IR_E_DEF, IR_R_DEF, IR_L_DEF] = stlp(100*GDPDefl(loc_start+1:loc_end-2),[pc(loc_start+1:loc_end-2,1:mpc)],Ztilde(1:end-2), ...
                  ProbRecession(loc_start:loc_end-1-2),lags,H,DTFP_UTIL(loc_start+1:loc_end-2));
            %             if iii == 3
            %                   IR_E = {-IR_E_T, -IR_E_G,-IR_E_C, -IR_E_U, -IR_E_W, -IR_E_H, -IR_E_CPI, -IR_E_INV};
            %                   IR_R = {-IR_R_T, -IR_R_G,-IR_R_C, -IR_R_U, -IR_R_W, -IR_R_H, -IR_R_CPI, -IR_R_INV};
            %                   IR_L = {-IR_L_T, -IR_L_G,-IR_L_C, -IR_L_U, -IR_L_W, -IR_L_H, -IR_L_CPI, -IR_L_INV};
            
            IR_E = {IR_E_T, IR_E_G,IR_E_C, IR_E_U, IR_E_W, IR_E_H, IR_E_CPI, IR_E_INV};
            IR_R = {IR_R_T, IR_R_G,IR_R_C, IR_R_U, IR_R_W, IR_R_H, IR_R_CPI, IR_R_INV};
            IR_L = {IR_L_T, IR_L_G,IR_L_C, IR_L_U, IR_L_W, IR_L_H, IR_L_CPI, IR_L_INV};
            
            nvar = length(varlist);
            n_row = 2;
            n_col = ceil(nvar/n_row);
            figure(3+iii)
            for j = 1: length(varlist)
                  s = subplot(n_row,n_col,j);
                  hold on
                  if j == 4 %|| j == 1 %|| j == 8 %do not take cumsum
                        q = plot([1:H]',IR_E{j}, '-r', 'linewidth', 3);
                        h = plot([1:H]',IR_R{j}, '--b','linewidth', 3);
                        l = plot([1:H]',IR_L{j}, '-ok','linewidth', 3);
                        plot([1:H]', 0*[1:H]', ':k');
                        set(gca,'TickLabelInterpreter','latex')
                        title(varlist{j},'interpreter', 'latex', 'fontsize', 12);
                  else
                        q = plot([1:H]',cumsum(IR_E{j}), '-r', 'linewidth', 3);
                        h = plot([1:H]',cumsum(IR_R{j}), '--b','linewidth', 3);
                        l = plot([1:H]',cumsum(IR_L{j}), '-ok','linewidth', 3);
                        plot([1:H]', 0*[1:H]', ':k');
                        set(gca,'TickLabelInterpreter','latex')
                        title(varlist{j},'interpreter', 'latex', 'fontsize', 14);
                  end
                  if j == 1
                        xlabel('Quarter','interpreter','latex','fontsize',12);
                        ylabel('\% deviation from s.s.','interpreter','latex','fontsize',12);
                  end
                  set(s, 'xlim', [1,H], 'ylim', ylim );
            end
            l=legend([q h l],{'Expansion','Recession','Linear'},'interpreter','latex');
            set(l, 'box','on', 'FontSize',13,'Orientation','horizontal','Position',[0.3 0.015 0.40 0.01]);
            
      end
      
      
end













