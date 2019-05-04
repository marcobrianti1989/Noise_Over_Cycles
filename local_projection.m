function [IRF,res,tuple_store,VDstore,degreeFreedom,nRegressor] = ...
      local_projection(y,x,u,lags,H,which_trend)

% H is the horizon of the irf
% y and x are (T,1); dependend and observable control variables
% u is exogenous regressor
% lags is how many lags of x and u we want to take in the regression
% regression: y(t+h) = alpha + B*u(t) + C1*u(t-1) + D1*y(t-1) + G1*x(t-1) + ...
y = detrend_func(y,which_trend);
for h = 1:H
      Y          = y(lags+h:end,:);
      X          = u(lags+1:end-h+1,:); %Plus 1 because when jj = 1 is contemporaneous
      if lags > 0 % which allows for controls
            for jj = 1:lags
                  if sum(sum(x.^2)) > 0 % Add into (big) X controls (small) x
                        X = [X, u(lags-jj+1:end-jj-h+1,:), y(lags-jj+1:end-jj-h+1,:), ...
                              x(lags-jj+1:end-jj-h+1,:)];
                        if h > 1 && jj == 1
                              X = [X, resh(2:end)];
                        end
                  else % no controls
                        X = [X, u(lags-jj+1:end-jj-h+1,:), y(lags-jj+1:end-jj-h+1,:)];
                        if h > 1 && jj == 1
                              X = [X, resh(2:end)];
                        end
                  end
            end
      end
      
      X                  = [X, ones(length(Y),1)];
      [T,n]              = size(X);
      degreeFreedom(h)   = T - n;
      nRegressor(h)      = n;
      tuple_store{h}     = [Y X];
      B                  = X'*X\(X'*Y);
      IRF(h)             = B(1);
      resh               = Y - X*B;
      res{h}             = resh;
      Rsquared(h)        = 1 - var(res{h})/var(Y);
      
      trunc = 25;
      
      XVD                       = X(:,2:end);
      BVD                       = XVD'*XVD\(XVD'*Y);
      resVD                     = Y - XVD*BVD;
      if h == 1
            CC = ones(length(resVD),1);
      end
      if h == 1
            XXVD                      = [CC u(lags+h:end,:)];
      elseif h < trunc && h > 1
            XXVD                      = [XXVD(1:end-1,:) u(lags+h:end,:)];
      else
            XXVD                      = [XXVD(1:end-1,1) XXVD(1:end-1,3:end) u(lags+h:end,:)];
      end
      [~,~,~,~,Rq]              = regress(resVD,XXVD);
      T                         = length(resVD);
      RVD(h)                    = Rq(1); %1 - (1-Rq(1))*((T-1)/(T-h-1)); adj Rsquared
           
      
      %Variance Decomposition
      %Den0 = var(resh) + Num; Den1 = var(resh - IRF(1)*X(1+1:end,1)) + Num; Den2 = var(resh - IRF(1)*X(1+2:end,1) - IRF(2)*X(1+1:end-1 ,1)) + Num;
%       if h == 1
%             NUM          = IRF.^2*var(u);
%             DEN          = var(resh) + NUM;
%       else
%             if h == 2
%                   Xstore       = u(lags+h:end,1);
%                   IRFrep       = repmat(IRF(1:end-1),[size(Xstore,1),1]);
%                   NUM          = sum(IRF.^2)*var(u);
%             elseif h > 2 && h <= trunc
%                   Xstore       = [u(lags+h:end,1) Xstore(1:end-1,:)];
%                   IRFrep       = repmat(IRF(1:end-1),[size(Xstore,1),1]);
%                   NUM          = sum(IRF.^2)*var(u);
%             else
%                   Xstore       = [u(lags+h:end,1) Xstore(1:end-1,1:end-1)];
%                   IRFrep       = repmat(IRF(1:trunc-1),[size(Xstore,1),1]);
%                   NUM          = sum(IRF(1:trunc-1).^2)*var(u);
%             end
%             DEN          = var(resh - sum(IRFrep.*Xstore,2)) + NUM;
%       end
%       VDstore(h)   = NUM/DEN;
      
      clear X
      
end

VDstore                   = RVD;