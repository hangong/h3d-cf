function [ei,model] = cf_3D_H(oi,ri,opt)
%CF_3D-H estimates a 3-D color homography color transfer model.
%
%   CF_3D-H(OI,RI) returns the colour transfered source
%   image OI according to the target image RI.
%
%   Options (opt.*):
%   * downsampling_res: specify the resolution of downsampled images
%     for faster color transfer approximation. [] for disabling.
%   * use_denoise: enable deniose by bilaterial filtering.
%   * use_curve: enable tone-mapping esimation by polynomial models.
%   * plot_RGB: an option for ploting RGBs (debug use).
%

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Gong, H., Finlayson, G.D., Fisher, R.B. and Fang, F., 2017. 3D color
%   homography model for photo-realistic color transfer re-coding. The
%   Visual Computer, pp.1-11.

% default parameters
if ~exist('opt','var'), opt = []; end
if ~isfield(opt,'downsampling_res'), opt.downsampling_res = []; end
if ~isfield(opt,'use_denoise'), opt.use_denoise = true; end
if ~isfield(opt,'use_curve'), opt.use_curve = true; end
if ~isfield(opt,'plot_RGB'), opt.plot_RGB = false; end

sz = size(oi);
f1 = reshape(oi,[],3);
f2 = reshape(ri,[],3);

% downsampling test
if ~isempty(opt.downsampling_res)
    % downsampling
    soi = imresize(oi,opt.downsampling_res);
    sri = imresize(ri,opt.downsampling_res);
    sf1 = reshape(soi,[],3)';
    sf2 = reshape(sri,[],3)';
else
    % use full res-images
    sf1 = f1'; sf2= f2';
end

% ploting RGB correspondences
if opt.plot_RGB
    figure;
    scatter3(sf1(1,:),sf1(2,:),sf1(3,:),40,sf1','filled','MarkerEdgeColor','r'); hold on;
    scatter3(sf2(1,:),sf2(2,:),sf2(3,:),40,sf2','filled','MarkerEdgeColor','k');
    for i = 1:size(sf1,2)
        plot3([sf1(1,i);sf2(1,i)],[sf1(2,i);sf2(2,i)],[sf1(3,i);sf2(3,i)],'b');
    end
    axis equal; xlim([0,1]);ylim([0,1]);zlim([0,1]);
end

% estimate 3D homography
P = [sf1;ones(1,size(sf1,2))];
Q = [sf2;ones(1,size(sf2,2))];
msk = min(P,[],1)>1/255 & min(Q,[],1)>1/255;
model.H = uea_H_from_x_als(P(:,msk),Q(:,msk),10);

% apply 3D homography
pe = model.H*[sf1;ones(1,size(sf1,2))];
pe = bsxfun(@rdivide,pe(1:3,:),pe(4,:));
pe = min(max(pe,0),1);

% brightness transfer
meanpe = mean(pe(:,msk),1)'; % transformed brightness
meante = mean(sf2(:,msk),1)'; % target brightness

if opt.use_curve
    % estimate brightness transfer
    model.pp = cvfit(meanpe,meante,'quad'); % b-b mapping
else
    % histogram matching
    model.pp = cvfit(meanpe,meante,'hist'); % b-b mapping
end

%pe = max(min(pe,1),0);

% re-apply to a higher res image
pe = model.H*[f1';ones(1,size(f1,1))];
pe = bsxfun(@rdivide,pe(1:3,:),pe(4,:));
pe = min(max(pe,0),1);
n = size(pe,2);
meanpe = mean(pe,1)'; % transformed brightness
if opt.use_curve
    meanf = model.pp(1+round(meanpe*255));
else
    meanf = model.pp(1+round(meanpe*255));
end
meanf = max(meanf,0);
nd = meanf./meanpe; % convert brightness change to shading
nd(meanpe<1/255) = 1; % discard dark pixel shadings
D = sparse(1:n,1:n,nd(:),n,n);

ei = reshape(pe',sz);
ImD = full(reshape(diag(D),sz(1:2)));

if opt.use_denoise % denoise the shading field
    grey = rgb2gray(oi);
    ImD = bFilter(ImD,grey,0,1,12);
end

% Debug
%{
imwrite(ImD/2,'s.jpg');
%}

ei = min(max(ei.*ImD,0),1);

end

function [H,err,d] = uea_H_from_x_als(P,Q,max_iter,tol)

% [H,rms,pa] = uea_H_from_x_als(H0,p1,p2,max_iter,tol)
%
% Compute H using alternating least square

% An initial estimate of
% H is required, which would usually be obtained using
% vgg_H_from_x_linear. It is not necessary to precondition the
% supplied points.
%
% The format of the xs is
% [x1 x2 x3 ... xn ; 
%  y1 y2 y3 ... yn ;
%  w1 w2 w3 ... wn]

if nargin<3, max_iter = 10; end
if nargin<4, tol = 1e-20; end

[Nch,Npx] = size(P);

% definition for Graham
fP = max(P,1e-6); fQ = max(Q,1e-6);
N = fP;

errs = Inf(max_iter+1,1); % error history

% solve the homography using ALS
n_it = 1; d_err = Inf;
while ( n_it-1<max_iter && d_err>tol)
    n_it = n_it+1; % increase number of iteration

    d = SolveD1(N,fQ);

    P_d = fP.*repmat(d,[Nch,1]);
    cv=P_d*P_d'; mma=mean(diag(cv));
    M = fQ*P_d'/(cv++eye(Nch).*mma./5000);
    N = M*fP;

    NDiff = (N.*repmat(d,[Nch,1])-Q).^2; % difference
    errs(n_it) = mean(mean(NDiff)); % mean square error
    d_err = errs(n_it-1) - errs(n_it); % 1 order error
end

H = M;
err = errs(n_it);

%plot(errs); hold on;
%fprintf('ALS %d: %f\n',n_it,errs(n_it));

%figure; imagesc(reshape(diag(D),4,6));
end

function d = SolveD1(p,q)
    nCh = size(p,1);
    d = (ones(1,nCh)*(p.*q))./(ones(1,nCh)*(p.*p));
end

