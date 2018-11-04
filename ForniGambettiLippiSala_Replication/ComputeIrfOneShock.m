function [C, sh] = ComputeIrfOneShock(B,GR,rsh)
N=size(B,1);
h=size(B,3);
if nargin == 2,
    sh = [];
    if size(GR,2)==0
        C = ones(N,0,h);
    else
        for k = 1:h
            C(:,:,k) = B(:,:,k)*GR;
        end
    end
else
    if size(GR,2) == 0
        C = ones(N,0,h);
        sh = ones(size(rsh,1),1,0);
    else
        for k = 1:h
            C(:,:,k) = B(:,:,k)*GR;
        end
        sh = rsh*GR;
    end
end
