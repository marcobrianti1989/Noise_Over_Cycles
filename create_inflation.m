function Xinfl = create_inflation(X,n)

% This function take X = logCPI which is (T,nvar) and 
% take differencence (in logs are growth rates) over n period of times.
% If t is a quarter and n is 4 then you get year on year inflation.

[T, nvar] = size(X);
for ii = 1:T-n
      Xinfl(ii) = X(ii+n,:) - X(ii,:);
end
Xinfl     = [NaN(n,nvar); Xinfl'];

end