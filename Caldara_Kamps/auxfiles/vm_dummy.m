% PURPOSE: Generates dummy observations for a Minnesota Prior
% -----------------------------------------------------
% USAGE: vm_dummy
% -----------------------------------------------------
% NOTE: requires to run vm_spec.m
% -----------------------------------------------------
 
% tau    : Overall tightness
% d      : Scaling down the variance for the coefficients of a distant lag
% w      : Number of observations used for obtaining the prior for the 
%          covariance matrix of error terms (usually fixed to 1). 
% lambda : Tuning parameter for coefficients for constant.
% mu     : Tuning parameter for the covariance between coefficients*/



nobs    = size(YY,1)-T0;  % number of observations


%--------------------
% Dummy Observations                                   
%--------------------

tau	   =   0.2; 
d	   =   2.5;
w	   =   0;
lambda =   0.5;	 
mu	   =   0.5;

if do_Figure == 5      % We lowered these hyper-parameters when variables are non-stationary
                       % to minimize the influence of the prior on the long-run dynamics of the
                       % model. No attempt was made to optimally select the hyper-parameters.
    lambda =   0.01;	 
    mu	   =   0.01;
end

YY0     =   YY(1:T0,:);  
ybar    =   mean(YY0)';      
sbar    =   std(YY0)'; 
premom  =   [ybar sbar];


%------------------------------------------
% Generate matrices with dummy observations
%------------------------------------------

hyp = [tau; d; w; lambda; mu];
[YYdum, XXdum, breakss] = varprior_h(n,p,nex_,hyp,premom);


%--------------------
% Actual observations
%--------------------

YYact = YY(T0+1:T0+nobs,:);
XXact = zeros(nobs,n*p);

i = 1;

while (i <= p)
    XXact(:,(i-1)*n+1:i*n) = YY(T0-(i-1):T0+nobs-i,:);
    i = i+1;
end

XXact = [XXact ones(nobs,1)];