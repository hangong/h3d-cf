function [ei,H] = cf_2D_H(oi,ri,opt)
%CF_2D-H estimates a 2-D color homography color transfer model.
%
%   CF_2D_H(OI,RI) returns the colour transfered source
%   image OI according to the target image RI.
%
%   Options (opt.*):
%   * downsampling_res: specify the resolution of downsampled images
%     for faster color transfer approximation. [] for disabling.
%   * use_denoise: enable deniose by bilaterial filtering.
%   * use_curve: enable tone-mapping esimation by polynomial models.
%

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Gong, H., Finlayson, G.D., Fisher, R.B.: Recoding color transfer as a
%   color homography. In: British Machine Vision Conference. BMVA (2016)

if ~exist('opt','var'), opt = []; end
if ~isfield(opt,'use_denoise'), opt.use_denoise = true; end
if ~isfield(opt,'use_curve'), opt.use_curve = true; end

if nargin<3, mapping = true; end

sz = size(oi);

f1 = reshape(oi,[],3);
f2 = reshape(ri,[],3);
sf1 = f1'; sf2= f2';

ssz = sz;
% compute homography
P = sf1; Q = sf2;
P = bsxfun(@rdivide, P, sum(P,1));
Q = bsxfun(@rdivide, Q, sum(Q,1));
C = [1,0,0;0,1,0;1,1,1];
H = H_from_x_als(P,Q,20);

H = C\(C*H);

pe = H*sf1; % Apply homography
pe = max(pe,0);

D = SolveD1(pe,sf2);

if ~opt.use_curve
    % apply shading
    D_new = min(D,10); % avoid shading artefact
else
    d = diag(D);
    meanpe = mean(pe,1)';

    %d = min(d,10);
    pp = csfit(meanpe,d,50);
    nd = ppval(pp,meanpe);
    if opt.use_curve
        D_new = SolveD2(nd,ssz);
    else
        nPx = size(nd,1);
        D_new = spdiags(nd,0,nPx,nPx);
    end
end

% fix D, debug use.
%{
ImD = full(reshape(diag(D),sz(1:2)));
ImD = bFilter(ImD);
imwrite(ImD/2,'s.jpg');
D_new = sparse(1:numel(ImD),1:numel(ImD),reshape(ImD,[],1)');
%}

pe = pe*D_new;
ei = reshape(pe',sz);

% compare with direct bilaterial filtering
%{
for i = 1:3
    ei(:,:,i) = bFilter(ei(:,:,i));
end
imwrite(ei,'s.jpg');
%}

end

function D = SolveD1(p,q)
    [nCh,nPx] = size(p);
    d = (ones(1,nCh)*(p.*q))./(ones(1,nCh)*(p.*p));
    D = spdiags(d',0,nPx,nPx);
end

function D = SolveD2(nd,sz)

    nPx = size(nd,1);

    A1 = speye(nPx);

    % compute D
    M1 = ShadingDiff(sz(1:2));
    lambda = 1./mean(diag(M1));
    lambda = 0.1*lambda;
    D = (A1+lambda*M1)\(nd);

    D = spdiags(D,0,nPx,nPx);
end

function M2 = ShadingDiff(lsz)
% minimise an edge image

    nel = prod(lsz);
    snel = prod(lsz-2);

    ind = zeros(lsz); ind(:) = 1:nel;
    cdx = ind(2:lsz(1)-1,2:lsz(2)-1); % centre
    tdx = ind(1:lsz(1)-2,2:lsz(2)-1); % top
    bdx = ind(3:lsz(1),2:lsz(2)-1);   % bottom
    ldx = ind(2:lsz(1)-1,1:lsz(2)-2); % left
    rdx = ind(2:lsz(1)-1,3:lsz(2));   % right

    % flatten index
    cdx = cdx(:); tdx = tdx(:); bdx = bdx(:); ldx = ldx(:); rdx = rdx(:);

    sM = sparse(cdx,cdx,-4*ones(1,snel),nel,nel) + ... % centre 
         sparse(cdx,tdx,ones(1,snel),nel,nel) + ... % top
         sparse(cdx,bdx,ones(1,snel),nel,nel) + ... % bottom
         sparse(cdx,ldx,ones(1,snel),nel,nel) + ... % left
         sparse(cdx,rdx,ones(1,snel),nel,nel); % right
    M2 = sM'*sM;
end


function [H,err,pD] = H_from_x_als(p1,p2,max_iter,tol)

if nargin<3, max_iter = 50; end
if nargin<4, tol = 1e-20; end

[Nch,~] = size(p1);

% definition for Graham
P = max(p1,1e-6);
Q = max(p2,1e-6);
N = P;

errs = Inf(max_iter+1,1); % error history

% solve the homography using ALS
n_it = 1; d_err = Inf;
while ( n_it-1<max_iter && d_err>tol )
    n_it = n_it+1; % increase number of iteration

    D = SolveD1(N,Q);

    P_d = P*D;
    cv = P_d*P_d'; mma = mean(diag(cv));
    M = Q*P_d'/(cv++eye(Nch).*mma./1000000);
    N = M*P;

    NDiff = (N*D-Q).^2; % difference
    errs(n_it) = mean(mean(NDiff)); % mean square error
    d_err = errs(n_it-1) - errs(n_it); % 1 order error
end

H = M;
err = errs(n_it);

end

