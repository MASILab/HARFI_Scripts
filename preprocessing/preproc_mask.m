

NIfTI_path = '/home-local/software/NIfTI_20140122';
addpath( NIfTI_path );

MNI_T1_template = '/nfs/share5/clineci/REX_TBI_VUMC/HARFI/HARFI_Scripts/preprocessing/MNI152_T1_2mm_brain.nii.gz';

wrk_dir =  '/nfs/share5/clineci/REX_TBI_VUMC/HARFI/data/preproc';
mask_nii = 'orig_target_seg.nii.gz';

subjects = dir(wrk_dir);

for i = 1:size(subjects,1)
    s = subjects(i);
    if size(s.name,2) > 3 % skip current dir & parent dir ('.','..')
        fprintf('*****************************************\n');
        fprintf('\t\t%s\n',s.name);
        fprintf('*****************************************\n');
        data_dir = [wrk_dir,'/',s.name,'/'];
        REG_dir = [data_dir, 'proc/Resting_State/REG/'];
        
        % register segmentation to MNI
        mask_path = [data_dir, mask_nii];
        out_fn = [data_dir,'mask2mni'];
        cmd = ['flirt -in ' mask_path, ...
               ' -ref ' MNI_T1_template, ...
               ' -applyxfm -init ' REG_dir,'reg_subT1_to_maMNI.mat' ...
               ' -interp trilinear -out ',out_fn,'.nii.gz'];
        display(cmd); system(cmd);
        
        % open vol and make all labels 1 (make brain mask)
        info = niftiinfo([out_fn,'.nii.gz']);
        mask = niftiread(info);
        mask(mask ~= 0) = 1;
        niftiwrite(mask,out_fn,info,'Compressed',true);
        
    end

end

% EPI to MNI  
    