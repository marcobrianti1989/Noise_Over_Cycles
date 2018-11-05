function B = FAVARCholBoot(Data,Z,variables,k,h,nrepli)
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
    B( :, :, :, j)  = FAVARCholImp(X_boot,Z_boot,variables,k, h);
end



