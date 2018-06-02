function [post_BETA_median, post_SIGMA_median, ...
      post_BETA_distribution, post_SIGMA_distribution] = ...
      bayesian_OLS_Gibbs_Sampling(prior_BETA,prior_VARbeta,prior_s,prior_v,...
      Y,X,N_simul,n)
%This function estimate a bayesian multivariate OLS using Gibbs sampling
%procedure. We allow dependent variable Y to be (nvar,T).

%%%%% Inputs: %%%%%
%prior_BETA: prior of the regression coefficient. (k,nvar)
%prior_VARbeta: prior of the variance of BETA. (k,k)
%prior_s: first argument of the Inverse-Wishart Distributin (nvar,nvar)
%prior_v: second argument of the Inverse-Wishart Distribution (1,1)
%Y: dependent variable vector. (T,nvar)
%X: independent variable matrix (T,k)
%N: number of simulations
%n: scalar usuful to derive posterior of v.

%%%%% Outputs: %%%%%
%post_BETA_distribution: posteriors of the regression coefficient. (k,nvar,N_simul)
%post_SIGMA_distribution: posteriors of the variance of the shocks to y. (nvar,nvar,N_simul)
%post_BETA_median: median of the posteriors of the regression coefficient. (k,nvar)
%post_SIGMA_median: median of posteriors of the variance of the shocks to y. (nvar,nvar)

%Technical parameters
[T, nvar]   = size(Y);
[T2, k]     = size(X);
if (T - T2)^2 > 10^(-14)
      error('X and Y do not have the same number of observations')
end

%Standard frequency approach OLS
B_hat = (X'*X)^(-1)*(X'*Y);
SSR = (Y-X*B_hat)'*(Y-X*B_hat);

%Standard procedure to obatin posterios in the Gibbs sampling
post_s = prior_s + SSR;
post_v = n*T + prior_v;

%Create Matrices for posterios parameters
post_BETA_distribution = zeros(k,nvar,N_simul);
post_SIGMA_distribution = zeros(1,nvar,N_simul);

%Gibbs Sampling Procedure
for i_simul = 1:N_simul
      %Step (1) Draw sigma from IW
      post_SIGMA_distribution(:,:,i_simul)  = iwishrnd(post_s,post_v);
      %Step (2) Draw BETA conditioning on sigma
      post_BETA_distribution(:,:,i_simul)   = (X'*X/post_SIGMA_distribution(:,:,i_simul) + ...
            prior_VARbeta^(-1))^(-1)*(X'*Y/post_SIGMA_distribution(:,:,i_simul) + ...
            prior_VARbeta^(-1)*prior_BETA)...
            + diag((X'*X/post_SIGMA_distribution(:,:,i_simul) + ...
            prior_VARbeta^(-1))^(-1)).*randn(k,nvar);
end

post_BETA_median   = median(post_BETA_distribution,3);
post_SIGMA_median  = median(post_SIGMA_distribution,3);


end