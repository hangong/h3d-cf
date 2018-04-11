function [ei,T] = cf_MK(I0,I1)
%CF_MK estimates a linear Monge-Kantorovitch color transfer model.
%
%   CF_POLY(I0,I1) returns the colour transfered source
%   image I0 according to the target image I1.
%

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.
%   Copyright F. Pitie 2007

%   References:
%   The linear Monge-Kantorovitch linear colour mapping for
%   example-based colour transfer. F. Pitié and A. Kokaram (2007) In 4th
%   IEE European Conference on Visual Media Production (CVMP'07). London,
%   November.

if (ndims(I0)~=3)
    error('pictures must have 3 dimensions');
end

X0 = reshape(I0, [], size(I0,3));
X1 = reshape(I1, [], size(I1,3));

A = cov(X0);
B = cov(X1);

T = MKL(A, B);

p1 = [X0,ones(size(X0,1),1)]; % Nx4

os1 = eye(4); os2 = eye(4);
os1(4,1:3) = -mean(X0); % 1x3
os2(4,1:3) = mean(X1); % 1x3

T(4,4) = 1;
T = os1*T*os2;
T = T(1:4,1:3);

XR = p1*T;

ei = reshape(XR, size(I0));

function [T] = MKL(A, B)
    N = size(A,1);
    [Ua,Da2] = eig(A); 
    Da2 = diag(Da2); 
    Da2(Da2<0) = 0;
    Da = diag(sqrt(Da2 + eps));
    C = Da*Ua'*B*Ua*Da;
    [Uc,Dc2] = eig(C); 
    Dc2 = diag(Dc2);
    Dc2(Dc2<0) = 0;
    Dc = diag(sqrt(Dc2 + eps));
    Da_inv = diag(1./(diag(Da)));
    T = Ua*Da_inv*Uc*Dc*Uc'*Da_inv*Ua';
end

end
