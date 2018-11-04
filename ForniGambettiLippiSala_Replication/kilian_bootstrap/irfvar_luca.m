% IRFVAR.M
% Lutz Kilian
% University of Michigan
% April 1997

function [IRF]=irfvar_luca(A,SIGMA,p,q,h)


J=[eye(q,q) zeros(q,q*(p-1))];
IRF=(J*A^0*J'*chol(SIGMA)');

for i=1:h
	IRF(:,:,i+1) = J*A^i*J'*chol(SIGMA)';
end;

