function C = polynomialmatricesproduct(A , B, nlags)
[A1, A2, A3] = size(A);
[B1, B2, B3] = size(B);
if nargin == 2
    nlags = A3 + B3 - 1;
end
C = zeros(A1, B2,A3 + B3 - 1 );
for a = 0: A3-1
    for b = 0: B3-1
        C(:,:,a + b + 1) = C(:,:,a + b + 1) + A(:,:,a + 1)*B(:,:,b + 1);
    end
end
C = C(:,:,1:nlags);