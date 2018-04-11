%   'test_main.m' is the main evaluation script for testing vairous color
%   transfer approximation methods. The error between the original color
%   transfer output and an approximated one is measured in PSNR an SSIM.
%   Also, we show a one-way anova post hoc test to analyse the similarity
%   between different results.

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Gong, H., Finlayson, G.D., Fisher, R.B. and Fang, F., 2017. 3D color
%   homography model for photo-realistic color transfer re-coding. The
%   Visual Computer, pp.1-11.

% configuration
in_path = 'in_pair/'; % enchancement output path
out_path = 'out_pair/'; % output path
mpath = './methods/'; % colour transfer method path

% define colour enhancement methods
ap.Name = {'MK','Poly','2D_H','3D_H'}; % name of the methods.
cf.Name = {'nguyen', 'pitie','pouli11','reinhard'}; % name of the methods.
Nap = numel(ap.Name); % number of approximation methods
Ncf = numel(cf.Name); % number of color transfer methods

% discover all images
Nf = 200;

if ~exist(out_path,'dir'), mkdir(out_path); end

t_psnr = zeros(Ncf,Nap);
t_ssim = zeros(Ncf,Nap);
a_psnr = cell(1,Nf); a_ssim = a_psnr;
for i = 1:Nf
    im_i = im2double(imread([in_path,num2str(i),'_s.jpg'])); % input image i
    im_j = im2double(imread([in_path,num2str(i),'_t.jpg'])); % input image j

    Ie = cell(Ncf,Nap);
    t_psnr(:) = NaN;
    for i_cf = 1:Ncf
        out_f = sprintf('%d_%s.jpg',i,cf.Name{i_cf});
        f_out = [in_path,out_f];
        if ~exist(f_out,'file')
            % apply a color transfer
            addpath([mpath,cf.Name{i_cf}]); % add colour correction method path
            cf_h = str2func(['cf_',cf.Name{i_cf}]); % function handle
            oIe = cf_h(im_i,im_j); % original colour transfer result
            oIe = max(min(abs(oIe),1),0);
            imwrite(oIe, f_out);
            rmpath([mpath,cf.Name{i_cf}]);
        end
        oIe = im2double(imread(f_out));

        for i_ap = 1:Nap
            ap_h = str2func(['cf_',ap.Name{i_ap}]); % function handle
            Ie{i_cf,i_ap} = ap_h(im_i,oIe);
            f_ap = sprintf('%s%d_%s_%s.jpg',out_path,i,cf.Name{i_cf},ap.Name{i_ap});
            imwrite(Ie{i_cf,i_ap}, f_ap);
            t_psnr(i_cf,i_ap) = psnr(Ie{i_cf,i_ap},oIe);
            t_ssim(i_cf,i_ap) = ssim(Ie{i_cf,i_ap},oIe);
        end
    end

    % save PSNR SSIM
    a_psnr{1,i} = t_psnr;
    a_ssim{1,i} = t_ssim;
end

% reshaping
a_psnr = reshape(cell2mat(a_psnr),[Ncf,Nap,Nf]);
a_ssim = reshape(cell2mat(a_ssim),[Ncf,Nap,Nf]);

% get mean
m_psnr = mean(a_psnr,3);
m_ssim = mean(a_ssim,3);

% construct tables
T_psnr = array2table(m_psnr','RowNames',ap.Name,'VariableNames',cf.Name);
T_ssim = array2table(m_ssim','RowNames',ap.Name,'VariableNames',cf.Name);

disp('psnr');
disp(T_psnr);
disp('ssim');
disp(T_ssim);

% compute the anova analysis
flat_psnr = reshape(permute(a_psnr,[1,3,2]),[],Nap);
flat_ssim = reshape(permute(a_ssim,[1,3,2]),[],Nap);

% perform one-way anova
[p_psnr,tbl_psnr,stats_psnr] = anova1(flat_psnr);
[p_ssim,tbl_ssim,stats_ssim] = anova1(flat_psnr);

% get a post hoc analysis on the group means
c_psnr = multcompare(stats_psnr);
c_ssim = multcompare(stats_ssim);

header = {'Method_A', 'Method_B', 'p_value'};
T_psnr_p = table(ap.Name(c_psnr(:,1))',ap.Name(c_psnr(:,2))', ...
           c_psnr(:,6),'VariableNames',header);
T_ssim_p = table(ap.Name(c_ssim(:,1))',ap.Name(c_ssim(:,2))', ...
           c_ssim(:,6),'VariableNames',header);

disp('psnr p-value');
disp(T_psnr_p);
disp('ssim p-value');
disp(T_ssim_p);
