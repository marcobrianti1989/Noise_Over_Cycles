function [dataset_boot] = data_boot_local_projection_wrong(B,res,nsimul,...
      which_correction,q,reg)
% Inputs:
% B: Beta coefficient (1,nreg)
% res = residuals from original regression
% nsimul = number of simulations, i.e. how many bootstrapped datasets to
% generate
% which_correction = 'none' or 'blocks'
% q = for block drawing, the length of each block
% reg is a vector of regressors

% Given regression Y = XB + E, I randomly reassign errors E to XB ...
% to build a new series Y' using a completely random extraction ...
% without replacement or assigning them in blocks of length q

T = size(res,1); % time periods
nreg = size(B,2);

switch which_correction
      case 'none'
            dataset_boot = zeros(T,nsimul);
            for i_repeat = 1:nsimul
                  for i_boot = 1:T
                        res_boot = res(randperm(T,1),:);
                        yhat = reg(i_boot,:)*B;
                        ystar = yhat + res_boot;
                        dataset_boot(i_boot,i_repeat) = ystar;
                  end
            end
      case 'blocks'
            dataset_boot = zeros(T,nsimul);            
            for i_repeat = 1:nsimul
                  i_boot = 1;
                  while i_boot <= T % cycles thru T in steps of q
                        draw = randperm(T,1); % draw one number from 1 to T
                        if draw >= T - q 
                              draw = draw - q - 1;
                        end
                        j = 1;
                        while j <= q && i_boot <= T % ... and for that index, go thru all shocks in that block
                              yhat = reg(i_boot,:)*B;
                              ystar = yhat + res(draw+j-1);
                              dataset_boot(i_boot,i_repeat) = ystar;
                              i_boot = i_boot + 1;
                              j = j + 1;
                        end
                  end
            end
      otherwise
            error('Correction needs to be either "none", or "blocks".')
            
end

end
