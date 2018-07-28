function [IRF_E, IRF_R, IRF_L , res_uncond, Rsquared, BL_store, XL_store] ...
      = stlp(y,x,u,fz,lags,H,TFP)
% H is the horizon of the irf
% y and x are (T,1);
% regression: y(t+h) = A0 + PI^E*(1-F(z(t-1)))*u(t) + PI^R*F(z(t-1))*u(t) + ...
%                         + B1*u(t-1) + C1*y(t-1) + D1*x(t-1) + C2*y(t-2) + ...

FZ              = repmat(fz,1,size(x,2));
%res_uncond      = zeros(length(Y),H);
for h = 1:H
      Y          = y(h+lags:end,:);
      X          = u(lags+1:end-h+1,:); %make sure fz is lagged in the argument of the function
      XL         = u(lags+1:end-h+1,:); %unconditional OLS, not controlling for smooth transition
      X          = [X, u(lags+1:end-h+1,:).*fz(lags+1:end-h+1,1)];
      if lags > 0 %which allows for controls
            for jj = 1:lags
                  if sum(sum(x.^2)) > 0 %Add into (big) X controls (small) x
                        X  = [X, u(lags-jj+1:end-jj-h+1,:), ...
                              u(lags-jj+1:end-jj-h+1,:).*fz(lags+1:end-h+1,1), ...
                              y(lags-jj+1:end-jj-h+1,:), ...
                              y(lags-jj+1:end-jj-h+1,:).*fz(lags+1:end-h+1,1), ...
                              x(lags-jj+1:end-jj-h+1,:), ...
                              x(lags-jj+1:end-jj-h+1,:).*FZ(lags+1:end-h+1,:)];
                        XL = [XL, u(lags-jj+1:end-jj-h+1,:), y(lags-jj+1:end-jj-h+1,:), ...
                              x(lags-jj+1:end-jj-h+1,:)];
                  else %no controls (small) x to add to (big) X
                        X  = [X, u(lags-jj+1:end-jj-h+1,:), ...
                              u(lags-jj+1:end-jj-h+1,:).*fz(lags+1:end-h+1,1), ...
                              y(lags-jj+1:end-jj-h+1,:), ...
                              y(lags-jj+1:end-jj-h+1,:).*fz(lags+1:end-h+1,1)];
                        XL = [XL, u(lags-jj+1:end-jj-h+1,:), y(lags-jj+1:end-jj-h+1,:)];
                  end
            end
      end
      % Add the linear trend - [1:1:length(Y)]'
      trend = [1:1:length(Y)]';
      if nargin > 6 %control for contemporaneous TFP
            X      = [X,  ones(length(Y),1), trend, TFP(h+lags:end,:)];
            XL     = [XL, ones(length(Y),1), trend, TFP(h+lags:end,:)];
      else
            X      = [X,  ones(length(Y),1), trend];
            XL     = [XL, ones(length(Y),1), trend];
      end
      B               = X'*X\(X'*Y);
      BL              = XL'*XL\(XL'*Y); %unconditional OLS, not controlling for smooth transition F(z)
      IRF_E(h)        = B(1); %here we are fixing the probability over the IRFs - need to relax
      IRF_R(h)        = B(1) + B(2); %this is correct, it can be proved mathematically in 2 steps
      IRF_L(h)        = BL(1); %unconditional OLS, not controlling for smooth transition F(z)
      res_uncond{h}   = Y - XL*BL;
      Rsquared(h)     = 1 - var(res_uncond{h})/var(Y);
      BL_store(:,h)   = BL;
      XL_store{h}     = XL;
end

close 


