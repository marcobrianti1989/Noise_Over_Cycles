%% Spectral density
% Compute and plot spectral density from MA.
% in the near future you may want to compute the spectrum, i.e. consider
% the multivariate case.
% TODO: accommodate multiple IRFs
function[spect, periodicity] = spectrum(IRF)
[r, c, d]  = size(IRF);
if c > r
    H = c;
    IRF = permute(IRF,[2 1 3]);
else
    H = r;
end
Sigma  = 1; % Shock variance
step   = .01;  % check relation with Canova frequencies
omega = 0 : step: pi;
for k = 1 : d
for x = 1: size(omega,2)
    for j = 1 : H
        one(:,j) = (IRF(j,:,k)*exp(-i*(j-1)*omega(1,x))); %CHECK (j-1)
        two(:,j) = (IRF(j,:,k)'*exp(i*(j-1)*omega(1,x))); %transpose is positive
    end
    sp(x,:)= (sum(one,2))*Sigma*(sum(two,2)); %CHECK
    %cross_spectrum1(x,1) = spectrum(2,1,x);
    %cross_spectrum1(x,2) = spectrum(1,2,x);%correggi: va bene perche quando ne fai il trasposto lo transforma in positivo!!%
end
cn = 1/(2*sum(sp,1));
spect(:,k) = sp.*cn;
end
periodicity = 2*pi./omega;




