function [] = fitHARFI(DIRcc,maskDIR)
% runs HARFI fitting on resting state fMRI data
% 
% Inputs
%   DIRcc : directory to correlation coefficients (contains directory to
%   nifti, directory to mask, R, discrete, g, cor_norm) (.mat)
%
% Outputs:
%   Will append SH_corr_even4, SH_corr_odd4, SH_corr_even2, and
%   SH_corr_odd2 to DIRcc .mat file
% 
% Dependencies:spherical_harmonics
%
% addpath('/Volumes/schillkg/MATLAB/spherical_harmonics/');
% DIRcc = '/Volumes/GRAID/FTI/HARFI/DetermineProcessing/sub3_r3_ccMatrix.mat'
% run_df_HARFI_fit(DIRcc)
%
% Written by Kurt Schilling - 02/2018 - kurt.g.schilling.1@vanderbilt.edu

load(DIRcc)
sz = size(cor_norm)
bvecs_v = g';
mask = load_untouch_nii(maskDIR); 
mask = single(mask.img);

SH_corr_even4 = zeros(sz(1),sz(2),sz(3),15);
SH_corr_odd4 = zeros(sz(1),sz(2),sz(3),25);
SH_corr_norm_even4 = zeros(sz(1),sz(2),sz(3),15);
SH_corr_norm_odd4 = zeros(sz(1),sz(2),sz(3),25);

for i = 1:sz(1)
    disp(i)
    for j = 1:sz(2)
        for k = 1:sz(3)
            if mask(i,j,k) > 0;
                r2 = squeeze(cor_norm(i,j,k,:));
                
                % Even 4th order SH coefficients
                SHr1 = SH_Fit2(r2,bvecs_v,4,0.006,'real'); % 15 by 1 double
                SH_corr_even4(i,j,k,:) = SHr1;
                % Odd 4th order SH coefficients 
                SHr2 = SH_Fit_all(r2,bvecs_v,4,0.006,'real',1); % 25 by 1 double
                SH_corr_odd4(i,j,k,:) = SHr2;
                
                % min max normalized r2
                r2_norm = (r2 - min(r2))/(max(r2) - min(r2));
                
                % Even 4th order SH coefficients
                SHr3 = SH_Fit2(r2_norm,bvecs_v,4,0.006,'real'); % 15 by 1 double
                SH_corr_norm_even4(i,j,k,:) = SHr3;
                % Odd 4th order SH coefficients
                SHr4 = SH_Fit_all(r2_norm,bvecs_v,4,0.006,'real',1); % 25 by 1 double
                SH_corr_norm_odd4(i,j,k,:) = SHr4;
      
            end
        end
    end
end

save(DIRcc,'SH_corr_even4','SH_corr_odd4','SH_corr_norm_odd4','SH_corr_norm_even4','-append');
