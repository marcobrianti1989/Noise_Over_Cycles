%----------------------------------------------------------------------
% Compute impulse responses and tax elasticities for the
% proxy SVAR approach identified using narrative tax shocks. 
% This code closely follows Merterns and Ravn (2013, AER)
%----------------------------------------------------------------------

function [EETA_PROXY_T, IRF_PROXY_T] =  proxy_svar_fiscal(Sigma,M,U,n,p,F,J,Hor)

    nIV = size(M,2);
    MM = [ones(length(M),1) M];
    Phib = MM\U;
    Phib = Phib(2:end,:);
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
    SS = chol(SigmaTSigmaTp);       
    ffactor1 = [b11iSig;b21iSig]*SS'; 
    EETA_PROXY_T = b12ib22;

    Omega1 = [ffactor1;zeros((p-1)*n,size(ffactor1,2))];

    Ltilde  = vm_irf(F,J,ffactor1,Hor,n,Omega1);
    IRF_PROXY_T = -Ltilde(:,3)/SS;
end
