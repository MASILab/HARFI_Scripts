function [] = HARFI_Diff(cc1,cc2,maskDIR,name)
% for every voxel
% for every direction
% calc difference in cor_norm
% difference is cc1-cc2
% save as CorDiff in "name"
% loop through all voxels, fit SH coefficients

% cc1 = '/Volumes/GRAID/HARFI/tHARFI/Dataset1/HARFI_PASSAGE_ALL.mat'
% cc2 = '/Volumes/GRAID/HARFI/tHARFI/Dataset1/HARFI_FIXATION.mat'
% maskDIR = '/Volumes/GRAID/HARFI/tHARFI/Dataset1/swravn3_mask.nii'
% name = 'PASSAGE-FIXATION.mat'

mask = load_untouch_nii(maskDIR); 
mask = single(mask.img);

load(cc1)
grad = g;
c1 = cor_norm;
clear cor_norm

load(cc2);
c2 = cor_norm;
clear cor_norm

sz = size(c1);


cor_diff = zeros(sz(1),sz(2),sz(3),length(grad));

for i = 1:sz(1)
    disp(i);
    for j = 1:sz(2)
        for k = 1:sz(3)
            if mask(i,j,k) > 0;
                for gg = 1:length(grad)
                    cor_diff(i,j,k,gg) = 2*c1(i,j,k,gg)-2*c2(i,j,k,gg);
                end
            end
        end
    end
end

cor_diff_relu = cor_diff;
cor_diff_relu(cor_diff_relu<0) = 0;

save(name,'cor_diff_relu','cor_diff');

%%
sz = size(cor_diff_relu)
bvecs_v = grad';

%SH_corr_even4 = zeros(sz(1),sz(2),sz(3),15);
%SH_corr_norm_even4 = zeros(sz(1),sz(2),sz(3),15);
SH_corr_odd4 = zeros(sz(1),sz(2),sz(3),25);
SH_corr_norm_odd4 = zeros(sz(1),sz(2),sz(3),25);

SH_corr_odd8 = zeros(sz(1),sz(2),sz(3),81);
SH_corr_norm_odd8 = zeros(sz(1),sz(2),sz(3),81);

SH_corr_odd12 = zeros(sz(1),sz(2),sz(3),169);
SH_corr_norm_odd12 = zeros(sz(1),sz(2),sz(3),169);


for i = 1:sz(1)
    disp(i)
    for j = 1:sz(2)
        for k = 1:sz(3)
            if mask(i,j,k) > 0;
                r2 = squeeze(cor_diff_relu(i,j,k,:));
                
                % Odd 4th order SH coefficients 
                SH = SH_Fit_all(r2,bvecs_v,4,0.006,'real',1); % 25 by 1 double
                SH_corr_odd4(i,j,k,:) = SH;
                %figure; showSH(SH,4,100,'real',1)
                
                % Odd 8th order SH coefficients 
                SH = SH_Fit_all(r2,bvecs_v,8,0.006,'real',1); % 81 by 1 double
                SH_corr_odd8(i,j,k,:) = SH;
                %figure; showSH(SH,8,100,'real',1)

                % Odd 12th order SH coefficients 
                SH = SH_Fit_all(r2,bvecs_v,12,0.006,'real',1); % 169 by 1 double
                SH_corr_odd12(i,j,k,:) = SH;
                %figure; showSH(SH,12,100,'real',1)
                
                % min max normalized r2
                r2_norm = (r2 - min(r2))/(max(r2) - min(r2));
                
                % Odd 4th order SH coefficients
                SH = SH_Fit_all(r2_norm,bvecs_v,4,0.006,'real',1); % 25 by 1 double
                SH_corr_norm_odd4(i,j,k,:) = SH;
                %figure; showSH(SH,4,100,'real',1)
                
                % Odd 8th order SH coefficients
                SH = SH_Fit_all(r2_norm,bvecs_v,8,0.006,'real',1); % 81 by 1 double
                SH_corr_norm_odd8(i,j,k,:) = SH;
                %figure; showSH(SH,8,100,'real',1)

                % Odd 12th order SH coefficients
                SH = SH_Fit_all(r2_norm,bvecs_v,12,0.006,'real',1); % 169 by 1 double
                SH_corr_norm_odd12(i,j,k,:) = SH;
                %figure; showSH(SH,12,100,'real',1)

            end
        end
    end
end

save(name,'SH_corr_odd4','SH_corr_norm_odd4','SH_corr_odd8','SH_corr_norm_odd8','SH_corr_odd12','SH_corr_norm_odd12','-append');











