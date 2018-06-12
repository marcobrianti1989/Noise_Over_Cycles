function pc = get_principal_components(data)
% Get principal components (see Stock & Watson 2002, eq 6)
% data is (T, nvar)

[T,n] = size(data);
for i=1:n
      data(:,i) = zscore(data(:,i));
end

% Get VC matrix
sigma = 1/T*data'*data;

[V,lamb] = eig(sigma); % diag(lamb) are the eigenvalues, V the eigenvectors
V = real(V);
lamb = real(lamb);
pc = data*V./n;

end
