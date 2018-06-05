function F21Ord1 = vm_mult_simple(F,J,Hor,eeta,U,T,n,p,nP,GDPRATIO)

U1 = U(:,1:nP);
U2 = U(:, nP+1:end);
u1hat = U1 - U2*eeta;
b11iSig1 = (u1hat'*u1hat)\u1hat'*U1;
b21iSig1 = (u1hat'*u1hat)\u1hat'*U2;
S1S1 = (u1hat'*u1hat)/(T-n*p-1);

% Compute two Cholesky orderings for T and G
s1 = sqrt(S1S1(1,1));
a  = S1S1(2,1)/s1;
s2 = sqrt(S1S1(2,2)-a^2);
SigmaOrd1 = [s1 0 ; a s2];
s2 = sqrt(S1S1(2,2));
b  = S1S1(1,2)/s2;
s1 = sqrt(S1S1(1,1)-b^2);
SigmaOrd2 = [s1 b; 0 s2];
ffactorOrd1 = [b11iSig1;b21iSig1']*SigmaOrd1;
ffactorOrd2 = [b11iSig1;b21iSig1']*SigmaOrd2;

OmegaTemp  = [ffactorOrd1;zeros((p-1)*n,size(ffactorOrd1,2))];
LtildeOrd1 = vm_irf(F,J,ffactorOrd1,Hor,n,OmegaTemp);

OmegaTemp  = [ffactorOrd2;zeros((p-1)*n,size(ffactorOrd2,2))];
LtildeOrd2 = vm_irf(F,J,ffactorOrd2,Hor,n,OmegaTemp);

F21Ord1(:,1) = -squeeze(LtildeOrd1(:,3,1))/SigmaOrd1(1,1)/GDPRATIO(2);
F21Ord1(:,2) =  squeeze(LtildeOrd2(:,3,2))/SigmaOrd2(2,2)/GDPRATIO(3);


