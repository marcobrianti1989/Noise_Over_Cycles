function [IRF_E, IRF_R, IRF_L ] = stlp(y,x,u,fz,lags,H,TFP)
% H is the length of the irf
% y and x are Tx1;
% the regression is  y(t+h) = A0 + PI^E(1-F(z(t-1)))u(t) + PI^RF(z(t-1))u(t) + B1u(t-1) + C1y(t-1) + D1x(t-1) + C2y(t-2) + ...

FZ = repmat(fz,1,size(x,2));
for h = 1:H
      Y = y(h+lags:end,:);
      X = u(lags+1:end-h+1,:); %make sure fz is lagged in the argument of the function
      XL = u(lags+1:end-h+1,:);
      X = [X, u(lags+1:end-h+1,:).*fz(lags+1:end-h+1,1)];
      if lags > 0 %which allows for controls
            for jj = 1:lags;
                  if x(1)^2 > 0;
                        X  = [X,u(lags-jj+1:end-jj-h+1,:), u(lags-jj+1:end-jj-h+1,:).*fz(lags+1:end-h+1,1), ...
                              y(lags-jj+1:end-jj-h+1,:), y(lags-jj+1:end-jj-h+1,:).*fz(lags+1:end-h+1,1),...
                              x(lags-jj+1:end-jj-h+1,:),x(lags-jj+1:end-jj-h+1,:).*FZ(lags+1:end-h+1,:)];
                        XL = [XL, u(lags-jj+1:end-jj-h+1,:), y(lags-jj+1:end-jj-h+1,:),...
                              x(lags-jj+1:end-jj-h+1,:)];
                  else
                        X  = [X,u(lags-jj+1:end-jj-h+1,:), u(lags-jj+1:end-jj-h+1,:).*fz(lags+1:end-h+1,1),y(lags-jj+1:end-jj-h+1,:), y(lags-jj+1:end-jj-h+1,:).*fz(lags+1:end-h+1,1)];
                        XL = [XL, u(lags-jj+1:end-jj-h+1,:), y(lags-jj+1:end-jj-h+1,:)];
                  end
            end
      end
      if nargin > 6
            X = [X , ones(length(Y),1), [1:1:length(Y)]',TFP(h+lags:end,:)];%,[1:1:length(Y)]' this add a linear trend
            XL = [XL , ones(length(Y),1),[1:1:length(Y)]',TFP(h+lags:end,:)];
      else
            X = [X , ones(length(Y),1), [1:1:length(Y)]'];%,[1:1:length(Y)]' this add a linear trend
            XL = [XL , ones(length(Y),1),[1:1:length(Y)]'];
      end
      B = X'*X\(X'*Y);
      BL = XL'*XL\(XL'*Y);
      IRF_E(h) = B(1); %here we are fixing the probability over the IRFs - need to relax
      IRF_R(h) = B(1) + B(2);
      IRF_L(h) = BL(1);
end




