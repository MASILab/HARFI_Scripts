function run_HARFI_preproc(wrk_dir, anat_img, fmri_img, NIfTI_path, SPM12_path)
% run_HARFI_preproc Run preprocessing pipeline for HARFI 
%   

    fprintf('Adding paths\n')
    addpath( NIfTI_path )
    addpath( SPM12_path )

    subjects = dir(wrk_dir);
    fprintf('Found %d subjects\n',size(subjects,1)-2)

    f = fopen('errorlog.txt','w');
    
    % get unzipped file names (strip off '.gz' if it's there)
    anat_nii = anat_img(1:end-3); 
    fmri_nii = fmri_img(1:end-3); 
    
    for i = 1:size(subjects,1)
        s = subjects(i);
        if size(s.name,2) > 3 % skip current dir & parent dir ('.','..')
            fprintf('*****************************************\n');
            fprintf('\t\t%s\n',s.name);
            fprintf('*****************************************\n');
            data_dir = [wrk_dir,'/',s.name,'/'];
            % check files
            if ~ ((exist([data_dir,anat_img],'file') && exist([data_dir,fmri_img],'file')) || (exist([data_dir,anat_nii],'file') && exist([data_dir,fmri_nii],'file')))
                fprintf(f,'Subject %s does not have a required file (T1 or fMRI)\n',s.name);
                continue
            end

            % Unzip images
            if strcmp(anat_img(end-2:end),'.gz') && (~ exist([data_dir,anat_nii],'file')) 
                gz_cmd = ['gunzip -f ', data_dir, anat_img];
                fprintf(gz_cmd)
                system(gz_cmd)
            end
            if strcmp(fmri_img(end-2:end),'.gz') && (~ exist([data_dir,fmri_nii],'file')) 
                gz_cmd = ['gunzip -f ', data_dir, fmri_img];
                fprintf(gz_cmd)
                system(gz_cmd)
            end


            % Pre-processes the T1 anatomical image and registers T1->MNI
            preproc_anat([data_dir,anat_nii]);

            
            % Pre-processes the fMRI image (drop first volumes, apply motion correction)
            preproc_fmri([data_dir,fmri_nii],{'drop','motco'},'',1);
            % get name of saved preprocessed fMRI vol
            fmri_name = fmri_nii(1:end-4);
            proc_fmri = sprintf('%s/proc/%s/%s_dr_mo.nii', data_dir,fmri_name,fmri_name);

            
            % Aligns fMRI to T1 and MNI
            preproc_alignment(proc_fmri,[data_dir,'T1_reg'],0);

            
            % Applies transformation matrices to all EPI volumes
            epi_REG_dir = sprintf('%s/proc/%s/REG/', data_dir, fmri_name);
            T1_REG_dir = [data_dir, 'T1_reg/'];
            align_epiToMNI(proc_fmri, epi_REG_dir, T1_REG_dir);
            mni_fmri = sprintf('%s/proc/%s/EPI_toMNI.nii', data_dir);
            
            
            % run nuisance regression on registered EPI vol
            preproc_fmri([data_dir,mni_fmri],{'nuis'},'',1);
        end
    end
end

