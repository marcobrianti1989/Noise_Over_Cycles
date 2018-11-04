function [NewsImp, SurpriseImp, Nsh, Ssh,impp] = FAVARNewsImp(X, Z, variables,k, h)


q = length(variables);
[impp,chi,sh] = FAVARCholImp(X, Z, variables, k,h);
Ssh = sh(:,1);
SurpriseImp = impp(:,1,:);
options = optimset('MaxFunEval',2000000,'MaxIter',1000000,'TolX',1e-15);

% [thetap, m, flag] = fminsearch(@LongRunEffect,sign(impp(variables(1),2:end,end))'.*ones(q-1,1),options,squeeze(impp(variables(1),:,end)));
% disp(flag)

[thetap, m, flag] = fminsearch(@LongRunEffectZ,sign(impp(variables(1),2:end,end))'.*ones(q-1,1),options,squeeze(impp(variables(1),:,[2:4 end])));
zetap=[0;thetap/norm(thetap)];
[NewsImp, Nsh] = ComputeIrfOneShock(impp,zetap, sh); % news



