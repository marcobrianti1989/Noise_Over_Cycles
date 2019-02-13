function [GxZ, Cx] = autocov_GenFunc(X,lags,Z)

[T, N] = size(X);
if N > 1
      error('X must be a unique process')
end
GxZ = 0;
ii = 1;
for i = -lags:lags
      Cx(ii,1)    = (X(1+lags:end-lags) - mean(X))'*(X(1+i+lags:end-lags+i)- mean(X))/(T-2);
      GxZ         = GxZ + Cx(ii,1)*Z^(-i);
      ii = ii + 1;
end