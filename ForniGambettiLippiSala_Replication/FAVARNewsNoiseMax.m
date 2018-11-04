function [irfs,irfsb,irf,irfb,sigma_a,sigma_aboot,bidielle,sh,ssh,wr,irfchol] = FAVARNewsNoiseMax(X,Z,codeData,k,h,nrepli,variables,ind1,ind,H)

%disp('FAVARNewsNoiseCholMax');

% ind1: position of the variable hit by "news shock" (potential, G,...)
% ind : position of the "signal" variable (G exp, or confidence....)

codeData2 = codeData;
if nargin == 9
H = h;
end
%__________________ Point estimates__________________

%[irf,~, shocks] = FAVARCholImp(X, Z, variables,k, h);
[N, S,Nsh_s,Ssh_s,irf] = FAVARNewsImp(X, Z, variables,k, H);
shocks(:,ind1) = Ssh_s;
shocks(:,ind) = Nsh_s;
irfchol=irf;
irf(:,2,:) = N;
irf(:,1,:) = S;

%irf(:,variables(ind1),:) = -irf(:,variables(ind1),:);
irf2 = irf;
    
if codeData(variables(ind1)) == 0
    for j = 1:size(irf,2)
        irf(variables(ind1),j,:) = filter([1 -1],1,squeeze(irf(variables(ind1),j,:)));
    end
end
cc = squeeze(irf(variables(ind1),ind,:));

%squeeze(irf(1,2,1:10))
%cc(1:10)

r = roots(flipud(cc(2:H)));
%abs(r)
wr = r(abs(r)<1);
%wr
bidielle = [ 0 1 zeros(1,h-1)];
if ~isempty(wr)
    for j = 1:length(wr)
        %                    MA          AR
        bidielle = filter([-wr(j) 1],[1 -conj(wr(j))],bidielle);
    end
end

%bidielle
bidielle = real(bidielle);

s = sum(cc(1:H))/sum(squeeze(irf(variables(ind1),ind1,(1:H))));
alphapoint = atan(s);
sigma_a = sin(alphapoint);% sigma_eps / sigma_s
sigma_e = cos(alphapoint);% sigma_v / sigma_s
disp([sigma_a sigma_e s])
%pnorm = 1;

if codeData(variables(ind1)) == 0
    irf(variables(ind1),:,:) = cumsum(irf(variables(ind1),:,:),3);
end

matrix = zeros(2,2,h+1);
matrix(1,1,:) = bidielle*sigma_e;
matrix(1,2,:) = -bidielle*sigma_a;
matrix(2,1,1) = sigma_a;
matrix(2,2,1) = sigma_e;

% compute structural irf
%size(irf(:,[ind1 ind],:))
%size(matrix)
%h+1
irfstruc = polynomialmatricesproduct(irf(:,[ind1 ind],:),matrix,h+1);
irfstruc(:,1,:) = irfstruc(:,1,:)*sign(s);

irfs = irfstruc;
irf = irf2;

%sh = [shocks(:,[1:ind-1 ind+1:end])*pnorm shocks(:,ind)];
sh = [shocks(:,ind1) shocks(:,ind)];
w = filter([ 0 1],1, flipud(sh(:,1)));
if ~isempty(wr)
    for j = 1:length(wr)
        w = filter([-wr(j) 1],[1 -conj(wr(j))],w) ;
    end
end
w=flipud(real(w));

ssh(:,1) = sigma_e*w + sigma_a*sh(:,2);
ssh(:,2) = -sigma_a*w + sigma_e*sh(:,2);

%__________________ Bootstrap replications __________________

irfb = FAVARCholBootMax(X,Z,variables,k,h,nrepli);
% irfb = FAVARCholBootAfterBoot(X,Z,variables,k,h,nrepli);

%
lo = 0;
for i=1:nrepli
    %i
    %irfb(:,ind1,:,i) = -irfb(:,ind1,:,i);
    
    irfb2(:,:,:,i) = irfb(:,:,:,i);
    if codeData(variables(ind1)) == 0
        for j = 1:size(irfb2,2)
            irfb(variables(ind1),j,:,i) = filter([1 -1],1,squeeze(irfb(variables(ind1),j,:,i)));
        end
    end
    
    cc = squeeze(irfb(variables(ind1),ind,:,i));
        
    r = roots(flipud(cc(2:H)));
    wrb = r(abs(r)<1);
    bidielleb = [ 0 1 zeros(1,h-1)];
    if ~isempty(wrb)
        for j = 1:length(wrb)
            bidielleb = filter([-wrb(j) 1],[1 -conj(wrb(j))],bidielleb) ;
        end
    end
    bidielleb = real(bidielleb);
    
    s = sum(cc(1:H))/sum(squeeze(irfb(variables(ind1),ind1,1:H,i)));
    alphapoint = atan(s);
    sigma_aboot(i) = sin(alphapoint);
    sigma_eboot(i) = cos(alphapoint);
    pnormb = 1;    
    
    if codeData(variables(ind1)) == 0
        irfb(variables(ind1),:,:,i) = cumsum(irfb(variables(ind1),:,:,i),3);
    end
        
    matrix = zeros(2,2,h+1);
    matrix(1,1,:) = bidielleb*sigma_eboot(i);
    matrix(1,2,:) = -bidielleb*sigma_aboot(i);
    matrix(2,1,1) = sigma_aboot(i);
    matrix(2,2,1) = sigma_eboot(i);
    
    % compute structural irf
    irfstruc(:,:,:,i) = polynomialmatricesproduct(irfb(:,[ind1 ind],:,i),matrix,h+1);
    irfstruc(:,1,:,i) = irfstruc(:,1,:,i)*sign(s);
    
end
irfsb =irfstruc;

