%% Spectral density
% Compute and plot spectral density from MA.
% in the near future you may want to compute the spectrum, i.e. consider
% the multivariate case.

H      = length(IRF);
Sigma  = 1; % normalization -- CHECK
step   = ?  % check relation with Canova frequencies 
lambda = 0 : step: pi;

for x = 1: size(lambda,2)
    for j = 1 : H
        uno(:,j) = (IRF(:,j)*exp(-i*(j-1)*lambda(1,x)));
        due(:,j) = (IRF(:,j)'*exp(i*(j-1)*lambda(1,x))); %transpose is positive
    end
    spectrum(:,x)= (sum(uno,3))*Sigma*(sum(due,3));
    cross_spectrum1(x,1) = spectrum(2,1,x);
    cross_spectrum1(x,2) = spectrum(1,2,x);%correggi: va bene perche quando ne fai il trasposto lo transforma in positivo!!%
end