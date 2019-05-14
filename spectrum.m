function[spect, periodicity] = spectrum(IRF)

% Spectral density
% Compute spectral density from truncated MA representation.
% in the near future you may want to compute the spectrum, i.e. consider
% the multivariate case.
% TODO: accommodate multiple IRFs


[r, h, nsimul]  = size(IRF);
if  r == 1
      H = h;
      IRF = permute(IRF,[2 1 3]);
else
      H = r;
end
Sigma  = 1; % Shock variance
step   = .01;  % check relation with Canova frequencies
omega = 0 : step: pi;
for k = 1 : nsimul                 % loop over bootstrap simulations
      for x = 1 : size(omega,2)    % loop over frequencies
            for j = 1 : H          % loop over MA(L)
                  one(:,j) = (IRF(j,:,k)*exp(-1i*(j-1)*omega(1,x))); %CHECK (j-1)
                  two(:,j) = (IRF(j,:,k)'*exp(1i*(j-1)*omega(1,x))); %transpose is positive
            end
            sp(x,:) = (sum(one,2))*Sigma*(sum(two,2)); %CHECK
            %cross_spectrum1(x,1) = spectrum(2,1,x);
            %cross_spectrum1(x,2) = spectrum(1,2,x);%correggi: va bene perche quando ne fai il trasposto lo transforma in positivo!!%
      end
      % Store spect
      spect(:,k) = sp;
      if k/50 == ceil(k/50)
            disp(['Spectrum: simulation number ',num2str(k)])
      end
end

periodicity = 2*pi./omega;


end

