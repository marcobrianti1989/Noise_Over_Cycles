function B = FAVARCholBootMax(Data,Z,variables,k,h,nrepli)
[T N] = size(Data);
r = size(Z, 2);
B = zeros( N, r, h + 1,nrepli);
[VarPa,C,X,u] = VarParameters(Z,k,1);
W = [ones(T,1) Z];
AA = (W'*W)\W'*Data;
chi = W*AA;

Idio = Data - chi;

for j=1:nrepli
    Z_boot = GenerateNewSeries(VarPa,C,X,u,k);
    W_boot = [ones(T-k,1) Z_boot];
    X_boot = W_boot*AA +Idio(k+1:end,:);
    %B( :, :, :, j)  = FAVARCholImp(X_boot,Z_boot,variables,k, h);
    [N, S,Nsh_s,Ssh_s,irf] = FAVARNewsImp(X_boot, Z_boot, variables,k, h);
    %shocks(:,ind1) = Ssh_s;
    %shocks(:,ind) = Nsh_s;
    irf(:,2,:) = N;
    irf(:,1,:) = S;
    B(:,:,:,j) = irf;
end



