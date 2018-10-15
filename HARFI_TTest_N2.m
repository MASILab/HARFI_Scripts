function [] = HARFI_TTest_N2(A1,A2,B1,B2,maskDIR,name)
% for every voxel
% for every direction
% calc p-value/t-value
% save as Pmat and Tmat in "name"
% test that A > B
% loop through all voxels, fit SH coefficients

mask = load_untouch_nii(maskDIR); 
mask = single(mask.img);

load(A1)
grad = g;
a1 = cor_norm;
clear cor_norm

load(A2);
a2 = cor_norm;
clear cor_norm

load(B1);
b1 = cor_norm;
clear cor_norm

load(B2);
b2 = cor_norm;
clear cor_norm

sz = size(a1);


Pmat = zeros(sz(1),sz(2),sz(3),length(grad));
Tmat = zeros(sz(1),sz(2),sz(3),length(grad));

for i = 1:sz(1)
    disp(i);
    for j = 1:sz(2)
        for k = 1:sz(3)
            if mask(i,j,k) > 0;
                for gg = 1:length(grad)
                    [h,p,ci,stats] = ttest2([a1(i,j,k,gg) a2(i,j,k,gg)],[b1(i,j,k,gg) b2(i,j,k,gg)],'tail','right');
                    Pmat(i,j,k,gg) = p;
                    Tmat(i,j,k,gg) = stats.tstat;
                end
            end
        end
    end
end

save(name,'Pmat','Tmat');

%%
sz = size(cor_diff_relu)
bvecs_v = grad';

SH_corr_odd8 = zeros(sz(1),sz(2),sz(3),81);



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



