function detrended_data = lin_timetrend(data,q)

% Simple OLS with Time Trend and Constant
T              = length(data);
Y              = data;
trend          = linspace(1,T,T)';
trendq         = trend.^2;
if q == 0
      X              = [ones(T,1) trend];
%       disp('Variables are linearly detrended')
%       fprintf('\n')
else
      X              = [ones(T,1) trend trendq];
%       disp('Variables are quadratically detrended')
%       fprintf('\n')
end
B              = (X'*X)^(-1)*X'*Y;
detrended_data = Y - X*B;


end