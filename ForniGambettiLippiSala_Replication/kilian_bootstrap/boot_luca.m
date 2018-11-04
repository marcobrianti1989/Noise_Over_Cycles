% Adapted from the code BOOT.M, by Lutz Kilian

function [CI,CI_chol] = boot_luca(A,U,y,V,p,h,nrep,codeData,ind,ind1)

[t,q]=size(y);
y = y';
Y = y(:,p:t);
for i = 1:p-1
    Y = [Y; y(:,p-i:t-i)];
end;

Ur = zeros(q*p,t-p);
Yr = zeros(q*p,t-p+1);
IRFrmat = zeros(nrep,q^2*(h+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  start of bootstrap simulation                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create nboot bootstrap replications of pseudo data

for j=1:nrep
    
    pos=fix(rand(1,1)*(t-p+1))+1;
    Yr(:,1)=Y(:,pos);
    
    index=fix(rand(1,t-p)*(t-p))+1;
    Ur(:,2:t-p+1)=U(:,index);
    
    for i=2:t-p+1
        Yr(:,i)= V + A*Yr(:,i-1)+Ur(:,i);
    end;
    
    yr=[Yr(1:q,:)];
    for i=2:p
        yr=[Yr((i-1)*q+1:i*q,1) yr];
    end;
    yr=yr';
    yr=detrend(yr,0);
    
    pr=p;
    
    [Ar,SIGMAr]=olsvarc(yr,pr);
    
    if ~ any(abs(eig(Ar))>=1)
        [Ar]=asybc(Ar,SIGMAr,t,pr);
    end;
    
    irf = irfvar_luca(Ar,SIGMAr(1:q,1:q),p,q,h);% Choleski + Kilian correction
    CI_chol(:,:,:,j) = irf;
    
    %for j = 1:size(irf,2)
    %    irf(variables(ind1),j,:) = filter([1 -1],1,squeeze(irf(variables(ind1),j,:)));
    %end
    
    %for jj = 1:size(irf,2)
    %    irf(1,jj,:) = filter([1 -1],1,squeeze(irf(1,jj,:)));
    %end
    %cc = squeeze(irf(1,2,:));
        
    if codeData(ind1) == 0
        for J = 1:size(irf,2)
            irf(ind1,J,:) = filter([1 -1],1,squeeze(irf(ind1,J,:)));
        end
    end
    cc = squeeze(irf(ind1,ind,:));
    
    r = roots(flipud(cc(2:end)));
    wr = r(abs(r)<1);
    bidielle = [ 0 1 zeros(1,h-1)];
    if ~isempty(wr)
        for jj = 1:length(wr)
            bidielle = filter([-wr(jj) 1],[1 -conj(wr(jj))],bidielle) ;
        end
    end
    bidielle = real(bidielle);
    
    s = sum(cc(1:h))/sum(squeeze(irf(ind1,ind1,(1:h))));
    alphapoint = atan(s);
    sigma_a = sin(alphapoint);
    sigma_e = cos(alphapoint);
    %pnorm = 1;
    
    %irf(1,:,:) = cumsum(irf(1,:,:),3);
    
    if codeData(ind1) == 0
        irf(ind1,:,:) = cumsum(irf(ind1,:,:),3);
    end
    
    matrix = zeros(2,2,h+1);
    matrix(1,1,:) = bidielle*sigma_e;
    matrix(1,2,:) = -bidielle*sigma_a;
    matrix(2,1,1) = sigma_a;
    matrix(2,2,1) = sigma_e;
        
    % compute structural irf
    
    irfstruc = polynomialmatricesproduct(irf(:,[ind1 ind],:),matrix,h+1);
    irfstruc(:,1,:) = irfstruc(:,1,:)*sign(s);
    
    %irfs = CumImp(irfstruc, codeData);
    %irf = CumImp(irf, codeData);

    shocks = U;
    sh = [shocks(:,1) shocks(:,2)];
    w = filter([ 0 1],1, flipud(sh(:,1)));
    if ~isempty(wr)
        for jj = 1:length(wr)
            w = filter([-wr(jj) 1],[1 -conj(wr(jj))],w) ;
        end
    end
    w=flipud(real(w));
    
    ssh(:,1) = sigma_e*w + sigma_a*sh(:,2);
    ssh(:,2) = -sigma_a*w + sigma_e*sh(:,2);
    
    %size(irfstruc)
    CI(:,:,:,j) = irfstruc;
end;

% Calculate 90 perccent interval endpoints
%CI=prctile(IRFrmat,[5 95]);

