function [F21Ord1, eetaTG] = vm_mult_bp(F,J,Hor,eeta,U,T,n,p,GDPRATIO)

%---------------------------
%Compute spending multiplier
%---------------------------
gres    = U(:,2);
macrores = U(:, 3:end);
gshock = gres - macrores*eeta(:,2);
L0gSig1 = (gshock'*gshock)\(gshock'*U);
S1S1 = (gshock'*gshock)/(T-n*p-1);
SS1 = chol(S1S1,'lower');
ffactorg = L0gSig1'*SS1;
OmegaTemp  = [ffactorg;zeros((p-1)*n,size(ffactorg,2))];
Ltildeg = vm_irf(F,J,ffactorg,Hor,n,OmegaTemp);
F21Ord1(:,1) = squeeze(Ltildeg(:,3,1))/SS1(1,1)/GDPRATIO(3);

%-----------------------
%Compute tax multiplier
%-----------------------

tres = U(:,1) - macrores*eeta(:,1);
eetaTG = (gshock'*gres)\gshock'*tres;
tshock = tres - eetaTG*gres;
L0tSig1 = (tshock'*tshock)\tshock'*U;
S1S1 = (tshock'*tshock)/(T-n*p-1);
SS1 = chol(S1S1,'lower');
ffactort = L0tSig1'*SS1;
OmegaTemp  = [ffactort;zeros((p-1)*n,size(ffactort,2))];
Ltildet = vm_irf(F,J,ffactort,Hor,n,OmegaTemp);
F21Ord1(:,2) = -squeeze(Ltildet(:,3,1))/SS1(1,1)/GDPRATIO(2);



