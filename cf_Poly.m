function [ei,M] = cf_Poly(oi,ri)
%CF_POLY estimates a polynomial color transfer model.
%
%   CF_POLY(OI,RI) returns the colour transfered source
%   image OI according to the target image RI.
%

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Ilie, A., Welch, G.: Ensuring color consistency across multiple cameras.
%   In: IEEE International Conference on Computer Vision, vol. 2,
%   pp. 1268–1275. (2005)

p1 = reshape(oi,[],3)';
p2 = reshape(ri,[],3)';

Npx = size(p1,2);

Pp1 = [p1;p1.^2;ones(1,Npx)];

% compute poly CC matrix
M = p2/Pp1;
%disp(M);

pe = M*Pp1; % chromaticity array

ei = reshape(pe',size(oi));

end
