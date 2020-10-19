function [] = align_epiToMNI(epi_fn_or_prefix, epi_REG_dir, T1_REG_dir)
% align epi to MNI data through pre-computed transforms.
% (epi->T1 via FSL's epi_reg, and T1->MNI via FSL as well)
% also resamples it to the MNI 2mm template.
%
% If a "dir/prefix" (without extension) is specified, will align
% all files having that prefix. Otherwise, will just align the
% specified file.
% 
% inputs:
% -------------------- %
% * epi_fn_or_prefix: fullpath to image (*.nii) or prefix
% (dir/name. will align dir/name*)
% * epi_REG_dir folder from which we can pull the transform output
% by fsl epi_reg (here, "oneVol_HC_am_EPI2T1.mat")
%      -  and make sure that T1_bet is there, too.
% * T1_REG_dir is where the Ants MNI transforms live.
%
% Sarah version 3/30/2020


MNI_T1_template = '/home-local/INN_final_project/preproc/MNI152_T1_2mm_brain.nii.gz';

[fp,fn,ext] = fileparts(epi_fn_or_prefix);
files_to_do = {};
if isempty(ext) % prefix
    prefix = [fp,'/',fn];
    display(['input is a prefix. Aligning all files: ',prefix,'*']);
    dd = dir([prefix,'*.nii']); 
    files_to_do = {};
    for jj=1:length(dd)
        files_to_do{jj} = [dd(jj).folder,'/',dd(jj).name];
    end
else
    if exist(epi_fn_or_prefix,'file')
        files_to_do{1} = epi_fn_or_prefix;
    end
end


display('*applying* alignment of EPI to subject T1 ');
display('----------------------------------------------');
epi_ref = [epi_REG_dir,'oneVol_HC_am'];
T1_bet_file = [epi_REG_dir,'T1_bet.nii']; assert(exist(T1_bet_file,'file')==2)
premat_file = [epi_ref,'_EPI2T1.mat']; assert(exist(premat_file,'file')==2)

T1_files = {};
for jj=1:length(files_to_do)
    [inpath,innm,inext] = fileparts(files_to_do{jj});
    out_file = [inpath,'/','EPI_toT1.nii.gz'];
    if ~exist(out_file,'file')
        cmd = ['applywarp  ', ...
               ' --ref=',T1_bet_file, ...
               ' --in=',files_to_do{jj}, ...
               ' --out=',out_file, ...
               ' --premat=',premat_file, ...
               ' --interp=trilinear'];
        display(cmd); system(cmd);
    else
       display(['already done: ',out_file]); 
    end
        % display(cmd); system(cmd);
        T1_files{jj} = out_file;
end


display('*applying* alignment of EPI to MNI ');
display('----------------------------------------------');
MNI_temp = [T1_REG_dir,'MNI152_T1_2mm_brain.nii.gz']; assert(exist(MNI_temp,'file')==2)
premat_file = [epi_REG_dir,'reg_subT1_to_maMNI.mat']; assert(exist(premat_file,'file')==2)

MNI_files = {};
for jj=1:length(T1_files)
    [inpath,innm,inext] = fileparts(T1_files{jj});
    out_file = [inpath,'/','EPI_toMNI.nii.gz'];
    if ~exist(out_file,'file')
        cmd = ['applywarp  ', ...
               ' --ref=',MNI_temp, ...
               ' --in=',T1_files{jj}, ...
               ' --out=',out_file, ...
               ' --premat=',premat_file, ...
               ' --interp=trilinear'];
        display(cmd); system(cmd);
    else
       display(['already done: ',out_file]); 
    end
        % display(cmd); system(cmd);
        MNI_files{jj} = out_file;
end

end

                             
