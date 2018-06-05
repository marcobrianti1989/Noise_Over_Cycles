function [eeta, Ltilde_fiscal] =  proxy_svar_nonfiscal(Sigma,U,n,p,F,J,Hor,M,M1,GDPRATIO,do_Figure)

LC =chol(Sigma,'lower');
A0 = (LC')\eye(size(LC,1));
if size(M,2) == 1
    nIV  = size(M,2);
    Phib = M\U;
else
    nIV  = size(M,2) + size(M1,2);
    Phib1 = M\U;
    Phib2 = M1\U(size(M,1)-size(M1,1)+1:end,:);
    Phib  = [Phib1(1,:); Phib2; Phib1(2,:)];
end


Phib11  = Phib(1:nIV,1:nIV);
Phib21  = Phib(1:nIV,nIV+1:n);
b21ib11 = (Phib11\Phib21)';
Sig11   = Sigma(1:nIV,1:nIV);
Sig21   = Sigma(nIV+1:n,1:nIV);
Sig22   = Sigma(nIV+1:n,nIV+1:n);
ZZp     = b21ib11*Sig11*b21ib11'-(Sig21*b21ib11'+b21ib11*Sig21')+Sig22;
b12b12p = (Sig21- b21ib11*Sig11)'*(ZZp\(Sig21- b21ib11*Sig11));
b11b11p = Sig11-b12b12p;
b22b22p = Sig22+b21ib11*(b12b12p-Sig11)*b21ib11';
b12ib22   = ((Sig21- b21ib11*Sig11)'+b12b12p*b21ib11')/(b22b22p');
b11iSig = eye(nIV)/(eye(nIV)-b12ib22*b21ib11);
b21iSig = b21ib11*b11iSig;
SigmaTSigmaTp =b11iSig\b11b11p/b11iSig';
SS = chol(SigmaTSigmaTp,'lower');       
b22iSig2 = (eye(n-nIV) - b21ib11*b12ib22)\eye(n-nIV);
b12iSig2 = b11iSig*b12ib22;
SigmaT2SigmaT2p =b22iSig2\b22b22p/b22iSig2';
SS2 = chol(SigmaT2SigmaT2p,'lower');
ffactor1 = [b11iSig;b21iSig]*SS; 
ffactor2 = [b12iSig2;b22iSig2]*SS2;
ffactor  = [ffactor1 ffactor2];
Q = LC\ffactor;
Qhat = Q(:,1:nIV);

%-------------------------
% Identify spending shock 
%-------------------------


Z = zeros(n-nIV-1,n);
if nIV == 1
    Z(1:2,2:3) = eye(2); % (under simple fiscal rule)
    Z(3,5) = 1;
else
    Z(1,5) = 1; %(under assumption that spending does not respond to contemporaneously taxes)
end

if do_Figure == 4;
    Z(end-1:end,6:7) = eye(2);
end

Rx = [Z*A0; Qhat'];
Nx=null(Rx);
Xhat=randn(n,1);
Qhatg =Nx*Nx'*Xhat/norm(Nx'*Xhat);

%-------------------------
% Identify tax shock 
%-------------------------

if nIV == 1
    Z = zeros(n-nIV-1,n); % (under simple fiscal rule)
    Z(1:2,2:3) = eye(2);
    Z(3,4) = 1;
    Rx = [Z*A0; Qhat'];
else
    Rx = [Qhat'; Qhatg']; % (uniquely identified by orthogonality to the other four shocks)
end

if do_Figure == 4; % Imposes Cholesky ordering of news variables
    if nIV == 1
        Z = zeros(n-nIV-1,n);
        Z(1:2,2:3) = eye(2);
        Z(end-1:end,6:7) = eye(2);
        Z(3,4) = 1;
        Rx = [Z*A0; Qhat'];
    else
        Z = zeros(n-nIV-2,n);
        Z(end-1:end,6:7) = eye(2);
        Rx = [Z*A0; Qhat'; Qhatg'];
    end
end


Nx=null(Rx);
Xhat=randn(n,1);
Qhatt =Nx*Nx'*Xhat/norm(Nx'*Xhat);

Q = [Qhat Qhatg Qhatt];
A0proxy = A0*Q;
ffactor  = LC*Q;

% Computing impulse responses

if  do_Figure == 6  % Stores IRFs for different ordering of gov spending and taxes
    Omega1 = [ffactor;zeros((p-1)*n,size(ffactor,2))];
    Ltilde = vm_irf(F,J,ffactor,Hor,n,Omega1);
    Ltilde_fiscal(:,1) =  Ltilde(:,1,nIV+2)*A0proxy(5,nIV+2)/GDPRATIO(3);
    Ltilde_fiscal(:,2) = -Ltilde(:,1,nIV+1)*A0proxy(4,nIV+1)/GDPRATIO(2);   
    eeta   = [-A0proxy(:,nIV+2)/A0proxy(5,nIV+2) -A0proxy(:,nIV+1)/A0proxy(4,nIV+1)];
elseif do_Figure == 7 % Stores IRFs for different scaling of multiplier
    Omega1 = [ffactor;zeros((p-1)*n,size(ffactor,2))];
    Ltilde = vm_irf(F,J,ffactor,Hor,n,Omega1);
    Ltilde_fiscal(:,1) =  Ltilde(:,1,nIV+1)/Ltilde(1,4,nIV+1)/GDPRATIO(3);
    Ltilde_fiscal(:,2) = -Ltilde(:,1,nIV+2)/Ltilde(1,5,nIV+2)/GDPRATIO(2);
    eeta   = [-A0proxy(:,nIV+1)/A0proxy(4,nIV+1) -A0proxy(:,nIV+2)/A0proxy(5,nIV+2)];
else
    Omega1 = [ffactor;zeros((p-1)*n,size(ffactor,2))];
    Ltilde = vm_irf(F,J,ffactor,Hor,n,Omega1);
    Ltilde_fiscal(:,1) =  Ltilde(:,1,nIV+1)*A0proxy(4,nIV+1)/GDPRATIO(3);
    Ltilde_fiscal(:,2) = -Ltilde(:,1,nIV+2)*A0proxy(5,nIV+2)/GDPRATIO(2);
    eeta   = [-A0proxy(:,nIV+1)/A0proxy(4,nIV+1) -A0proxy(:,nIV+2)/A0proxy(5,nIV+2)];
end

end
