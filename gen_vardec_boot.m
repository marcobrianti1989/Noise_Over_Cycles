function [vardec, ub1_vardec, lb1_vardec, ub2_vardec, lb2_vardec] = ...
      gen_vardec_boot(gamma,gamma_boot,chol,chol_boot,B,B_boot,H,sig1,sig2)

% Tools
nvar          = size(gamma,1);
nshocks       = nvar;
nsimul        = size(B_boot,3);
perc_up1      = ceil(nsimul*sig1); % the upper percentile of bootstrapped responses for CI
perc_low1     = floor(nsimul*(1-sig1)); % the lower percentile of bootstrapped responses for CI
perc_up2      = ceil(nsimul*sig2); % same for 2nd CI
perc_low2     = floor(nsimul*(1-sig2)); % same for 2nd CI
m             = linspace(1,H,H);
% Get variance Decomposition
if size(gamma,1) == size(gamma,2)
      impact_vardec = chol;
else
      N                     = null(gamma');
      D_null                = [gamma N];
      impact_vardec         = chol*D_null;
end
[IRF_vardec,~,~,~,~]  = genIRFs(impact_vardec,0,B,0,H,sig1,sig2);
for im = 1:length(m)
      vardec(:,im,:)     = gen_vardecomp(IRF_vardec,m(im),H);
end
% Get bootstrapped vardec nsimul times
for i_simul=1:nsimul
      if size(gamma,1) == size(gamma,2)
            impact_vardec_boot(:,:,i_simul) = chol_boot(:,:,i_simul);
      else
            N_boot(:,:,i_simul)             = null(gamma_boot(:,:,i_simul)');
            D_null_boot(:,:,i_simul)        = [gamma_boot(:,:,i_simul) N_boot(:,:,i_simul)];
            impact_vardec_boot(:,:,i_simul) = chol_boot(:,:,i_simul)*D_null_boot(:,:,i_simul); % where A is the chol.
      end
      [IRF_vardec_boot(:,:,:,i_simul),~,~,~,~] = ...
            genIRFs(impact_vardec_boot(:,:,i_simul),0,B_boot(:,:,i_simul),0,H,sig1,sig2);
      for im = 1:length(m)
            vardec_boot(:,im,:,i_simul) = gen_vardecomp(IRF_vardec_boot(:,:,:,i_simul),m(im),H);
      end
end

% Sort bootstrap IRFs and set lower and upper bounds
for i_shocks = 1:nshocks
      vardec_boot(:,:,i_shocks,:) = sort(vardec_boot(:,:,i_shocks,:),4);
      ub1_vardec(:,:,i_shocks)    = vardec_boot(:,:,i_shocks,perc_up1);
      lb1_vardec(:,:,i_shocks)    = vardec_boot(:,:,i_shocks,perc_low1);
      ub2_vardec(:,:,i_shocks)    = vardec_boot(:,:,i_shocks,perc_up2);
      lb2_vardec(:,:,i_shocks)    = vardec_boot(:,:,i_shocks,perc_low2);
end


end