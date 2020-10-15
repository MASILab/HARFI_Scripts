function [] = calcHARFI(DIRnii, DIRmask, R, discrete, DIRout, tr, N)
% runs HARFI on resting state fMRI data
% 
% Inputs
%   DIRnii: directory to resting state image (.nii)
%   DIRmask: directory to mask image (.nii)
%   R: radius of integration (voxels)
%   discrete: discrete steps in integration (100)
%   DIRout: directory for output (.mat)
%   tr: TR (s)
%   N: number of volumes to remove from beginning (recommend 20)
%
% Outputs
%   DIRout: .mat file with cor_norm, R, discrete, DIRnii, DIRmask, g
% 
% Dependencies: NIFTI package, spm12, YuruiCode (maybe add to end),
% Jones60.mat
%
% Written by Kurt Schilling - 01/2018 - kurt.g.schilling.1@vanderbilt.edu

%addpath('/Volumes/schillkg/MATLAB/spm12/');
%addpath('/Volumes/schillkg/MATLAB/NIFTI_20130306/');

% open FMRI
nii = load_untouch_nii(DIRnii); fmri = single(nii.img);
sz = size(fmri); aa = sz(1); bb = sz(2); cc = sz(3); dd= sz(4); nvols = sz(4);

% open mask
nii2 = load_untouch_nii(DIRmask); mask = single(nii2.img);

% mask
% for ii=1:nvols
%     fmri(:,:,:,ii) = squeeze(fmri(:,:,:,ii)).*mask;
% end

% remove first N and last M volumns
%N = 20; 
M = 0;
nvols = nvols - N - M;
fmri(:,:,:,1:N) = []; fmri(:,:,:,end-M+1:end) = [];
timepts = ((1:nvols)+N)*tr;

% plot global signal
sum4d = squeeze(sum(sum(sum(fmri,1),2),3));
figure(3); plot(timepts, sum4d(1:nvols),'b'); xlabel('time [s]');
%print('-dpng', '-r100', 'globalts.png');

% Detrend
disp('Detrending...')
gmwm_msk = mask; gmwm_ind = find(gmwm_msk>0);
[rows, cols, vols] = ind2sub(sz(1:3), gmwm_ind);
for nn = 1:length(rows) % per valid voxel
    row = rows(nn); col = cols(nn); vol = vols(nn); % voxel coord
    ts_vox = squeeze(fmri(row,col,vol,:));
    % fmri(row,col,vol,:) = detrendnonlin(ts_vox); % nonlinear detrend
    fmri(row,col,vol,:) = detrend(ts_vox); % linear detrend
end

% plot global signal
sum4d = squeeze(sum(sum(sum(fmri,1),2),3));
figure(4); plot(timepts, sum4d(1:nvols),'b'); xlabel('time [s]');
%print('-dpng', '-r100', 'globalts_detrend.png');

% Temperal filter 
disp('Temporal Filtering...')
fmri = reshape(fmri, [sz(1) sz(2) sz(3) nvols]);
[bb1, aa1] = mylpfilt2_ini(tr, 0.015, 0.01, 3, 60); % design filter, high pass
[bb2, aa2] = mylpfilt2_ini(tr, 0.1, 0.11, 3, 60); % deisng filter, low passfreqz(bb,aa,4096,1/tr); 
for nn = 1:length(rows) % per valid voxel
    row = rows(nn); col = cols(nn); vol = vols(nn); % voxel coord
    % load detrended signal
    ts_vox = double(squeeze(fmri(row,col,vol,:))); % single data type, Nx1
    
    if nn==1000 % display detrended data
        display(['1000th voxel- row:' num2str(row) ' col:' num2str(col) ' vol:' num2str(vol)]);
        figure(4); subplot(3,2,1); plot(timepts, ts_vox); xlim([timepts(1) timepts(end)]); title('after detrend (1000th voxel)');
        fft_vox = abs(fft(ts_vox));
        subplot(3,2,2); plot(1/tr*(0:(nvols/2))/nvols, fft_vox(1:(nvols/2+1))); xlabel('freqency (Hz)');
    end
    
    % high pass filtering
    nzeros = 20;
    ts_vox = cat(1, zeros([nzeros 1]), ts_vox, zeros([nzeros 1])); % zero padding
    ts_vox = filtfilt(bb1, aa1, ts_vox); % highpass filtering
    ts_vox(1:nzeros) = []; ts_vox(end-nzeros+1:end) = [];
    
    if nn==1000 % display high pass filtered data
        subplot(3,2,3); plot(timepts, ts_vox); xlim([timepts(1) timepts(end)]); title('after high pass filtering');
        fft_vox = abs(fft(ts_vox));
        subplot(3,2,4); plot(1/tr*(0:(nvols/2))/nvols, fft_vox(1:(nvols/2+1))); xlabel('freqency (Hz)');
    end
    
    % low pass filtering
    ts_vox = cat(1, zeros([nzeros 1]), ts_vox, zeros([nzeros 1])); % zero padding
    ts_vox = filtfilt(bb2, aa2, ts_vox); % lowpass filtering
    ts_vox(1:nzeros) = []; ts_vox(end-nzeros+1:end) = [];
    
    if nn==1000 % display low pass filtered data
        subplot(3,2,5); plot(timepts, ts_vox); xlim([timepts(1) timepts(end)]); title('after low pass filtering');
        fft_vox = abs(fft(ts_vox));
        subplot(3,2,6); plot(1/tr*(0:(nvols/2))/nvols, fft_vox(1:(nvols/2+1))); xlabel('freqency (Hz)');
        %print('-dpng', '-r100', 'smhplp.png');
    end
    
    fmri(row,col,vol,:) = single(ts_vox);
end

sz = size(fmri); aa = sz(1); bb = sz(2); cc = sz(3); dd= sz(4); 

% open direction schemes
load('Jones60.mat','g');
cor_norm = zeros(sz(1),sz(2),sz(3),length(g));

% integration distances
integrand = linspace(0,R,discrete);

x_m = integrand'*g(1,:); % matrix size = discrete X num of directions
y_m = integrand'*g(2,:);
z_m = integrand'*g(3,:);

disp('Calculating Correlations...')
for i = R-1:sz(1)-(R-1);
    disp(i)
    for j =  R-1:sz(2)-(R-1);
        for k =  R-1:sz(3)-(R-1);
            if mask(i,j,k) > 0;
                A = squeeze(fmri(i,j,k,:));
                B_m = zeros([discrete length(g) dd]);
                parfor kk = 1:dd % for each time point
                    B_m(:,:,kk) = interpn(1:aa,1:bb,1:cc,squeeze(fmri(:,:,:,kk)), i+x_m, j+y_m, k+z_m); % 3D matrix   
                end % parfor
                B_m = reshape(B_m, [discrete*length(g) dd]); % 3D -> 2D: each row is the fMRI fluctuation at one interplated position
                RR_m = corr(A, B_m'); % correlate vector A with each row vector of B_m; 
                RR_m = reshape(RR_m, [discrete length(g)]);
                cor_norm(i,j,k,:) = sum(abs(RR_m))/discrete; % sum of each column of RR_m;
            end % mask
        end % k
    end % j
end % i

save(DIRout,'cor_norm', 'R', 'discrete', 'DIRnii', 'DIRmask', 'g');


function [b, a] = mylpfilt2_ini(TR, Wp, Ws, Rp, Rs)

fs = 1/TR; % fs/2 = Nyquist freq
Wp = Wp/(fs/2); % Passband corner freq: normalized to Nyquist freq
Ws = Ws/(fs/2); % Stopband corner freq: normalized to Nyquist freq
[n, Wn] = cheb2ord(Wp, Ws, Rp, Rs);
if length(Wp)==1 && (Wp > Ws) % high pass
    [b, a] = cheby2(n, Rs, Wn, 'high');
else % low pass or band pass
    [b, a] = cheby2(n, Rs, Wn);
end
