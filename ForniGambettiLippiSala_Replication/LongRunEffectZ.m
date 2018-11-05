function lr = LongRunEffectZ(theta,C);
% zeta=[0;theta/norm(theta)];
% lr=-(C*zeta);
zeta=[0;theta/norm(theta)];
lr=-(C(:,end)'*zeta);
