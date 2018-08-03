function [Yboot Xboot] = bootstrap_LP(Y,X,nsimul,CI)

% Inputs:
% Y which is the projected variable at horizon h. Length is T
% X is the matrix of regressor. Length is T
% nsimul = number of simulations, i.e. how many bootstrapped datasets to generate
% As block_size we use approximately half of the dataset
% CI is confidence interval. Say 95% for example

% Given all the possible combinations of Y and X we extract only a subset
% of consecutive Y and X of length block_size

T                  = length(Y);
l                  = floor(T*0.3);
for isimul         = 1:nsimul
      draw               = randi(T-l,1);
      Yboot(:,isimul)    = Y(draw:draw+l);
      Xboot(:,:,isimul)  = X(draw:draw+l,:);
end








end