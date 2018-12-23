function [IRF,res,tuple_store,VDstore] = local_projection(y,x,u,lags,H,which_trend)

% H is the horizon of the irf
% y and x are (T,1); dependend and observable control variables
% u is exogenous regressor
% lags is how many lags of x and u we want to take in the regression
% regression: y(t+h) = alpha + B*u(t) + C1*u(t-1) + D1*y(t-1) + G1*x(t-1) + ...
NUM = 0;
for h = 1:H
      Y          = y(lags+h:end,:);
      X          = u(lags+1:end-h+1,:); %Plus 1 because when jj = 1 is contemporaneous
      if lags > 0 % which allows for controls
            for jj = 1:lags
                  if sum(sum(x.^2)) > 0 % Add into (big) X controls (small) x
                        X = [X, u(lags-jj+1:end-jj-h+1,:), y(lags-jj+1:end-jj-h+1,:), ...
                              x(lags-jj+1:end-jj-h+1,:)];
                        if h > 1 && jj == 1
                              X = [X, resh(1:end-1)];
                        end
                  else % no controls
                        X = [X, u(lags-jj+1:end-jj-h+1,:), y(lags-jj+1:end-jj-h+1,:)];
                        if h > 1 && jj == 1
                              X = [X, resh(1:end-1)];
                        end
                  end
            end
      end
      switch which_trend
            case 'BPfilter'
                  Y     = bpass(Y,4,32);
            case 'HPfilter'
                  [~,Y] = hpfilter(Y,1600);
            case 'linear'
                  trend = [1:1:length(Y)]';
            case 'quadratic'
                  trend = [1:1:length(Y)]';
                  trend = [trend trend.^2];
                  X     = [X trend];
            case 'none'
                  Y     = Y;
      end
      X                  = [X, ones(length(Y),1)];
      tuple_store{h}     = [Y X];
      B                  = X'*X\(X'*Y);
      IRF(h)             = B(1);
      resh               = Y - X*B;
      res{h}             = resh;
      Rsquared(h)        = 1 - var(res{h})/var(Y);
      
      % Variance Decomposition
      NUM                = NUM + B(1)^2;
      %Den0 = var(resh) + Num; Den1 = var(resh - IRF(1)*X(1+1:end,1)) + Num; Den2 = var(resh - IRF(1)*X(1+2:end,1) - IRF(2)*X(1+1:end-1 ,1)) + Num;
      if h == 1
            Xstore       = u(lags+1:end-h+1,1);
            DEN          = var(resh) + NUM;
      else
            Xstore       = [u(lags+h:end,1) Xstore(1:end-1,:)];
            IRFrep       = repmat(IRF,[size(Xstore,1),1]);
            DEN          = var(resh - sum(IRFrep.*Xstore,2)) + NUM;
      end
      VDstore(h)   = NUM/DEN;
end