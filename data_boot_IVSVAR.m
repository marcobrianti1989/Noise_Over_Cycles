function [dataset_boot, Z_boot] = data_boot_IVSVAR(B,nburn,res,Z,nsimul,which_correction)
% Inputs:
% B: Beta coefficient matrix from VAR (nvar*nlags + 1 , nvar)
% nburn = number of burn-in periods
% nlag = number of lags in VAR
% res = residuals from original VAR (Txnvar)
% nsimul = number of simulations, i.e. how many bootstrapped datasets to
% generate
% q = for block drawing, the length of each block

T  = size(res,1); % time periods
nvar = size(B,2);
if size(B,1)/size(B,2) == floor(size(B,1)/size(B,2))
      nlag = (size(B,1))/nvar;
      if ceil(nlag) == nlag
      else
            error('nlag number is not an integer')
      end
      B = [zeros(1,nvar); B];
else
      nlag = (size(B,1)-1)/nvar;
      if ceil(nlag) == nlag
      else
            error('nlag number is not an integer')
      end
end
const = 1;

switch which_correction
      case 'none'
            dataset_boot = zeros(T+nburn,nvar,nsimul);            
            for i_repeat = 1:nsimul
                  reg_new = zeros(1,nlag*nvar);
                  for i_boot = 1:T+nburn
                        flag                            = randperm(T,1);
                        res_boot                        = res(flag,:);
                        yhat                            = [const, reg_new]*B;
                        ystar                           = yhat + res_boot;                        
                        reg_new                         = [ystar reg_new(1:end-nvar)];
                        dataset_boot(i_boot,:,i_repeat) = ystar;
                        Z_boot(i_boot,:,i_repeat)       = Z(flag,:);                        
                  end
            end
      otherwise
            error('No Blocks correction so far.')
end

dataset_boot = dataset_boot(nburn+1:end,:,:);
Z_boot       = Z_boot(nburn+1:end,:);

end
