function run_HARFI_preproc(wrk_dir, anat_img, fmri_img, MNI_path, NIfTI_path, SPM12_path)
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
    fmri_name = fmri_nii(1:end-4);
    
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
            T1_REG_dir = [data_dir, 'T1_reg/']; % output dir
            anat_bet = [T1_REG_dir 'T1_bet.nii']; % output vol
            if ~exist(anat_bet,'file')
                preproc_anat([data_dir,anat_nii], MNI_path);
            end
 
            
            % Pre-processes the fMRI image (drop first volumes, apply motion correction)
            proc_fmri = sprintf('%sproc/%s/%s_dr_mo.nii', data_dir,fmri_name,fmri_name); % output vol
            motpar_fn = sprintf('%sproc/%s/%s_dr.volreg_par', data_dir,fmri_name,fmri_name); % motion param file name
            if ~( exist(proc_fmri,'file') && exist(motpar_fn,'file') )
                preproc_fmri([data_dir,fmri_nii],{'drop','motco'},'',1);
            end
            
            
            % Aligns fMRI to T1 and MNI
            ref_fmri = sprintf('%sproc/%s/REG/reg_toMNI_EPI.nii.gz', data_dir,fmri_name); % output vol
            if ~exist(ref_fmri,'file')
                preproc_alignment(proc_fmri,T1_REG_dir,0,MNI_path);
            end

            
            % Applies transformation matrices to all EPI volumes
            mni_fmri = sprintf('%sproc/%s/EPI_toMNI.nii', data_dir, fmri_name); % output vol
            if ~exist(mni_fmri,'file')
                epi_REG_dir = sprintf('%sproc/%s/REG/', data_dir, fmri_name);
                align_epiToMNI(proc_fmri, epi_REG_dir, T1_REG_dir);
                % Unzip registered EPI
                gz_cmd = ['gunzip -f ', mni_fmri,'.gz'];
                fprintf(gz_cmd)
                system(gz_cmd)
            end
           
            
            
            % Run nuisance regression on registered EPI vol
            fin_fmri = sprintf('%sproc/%s/EPI_toMNI_nr.nii', data_dir, fmri_name); % output vol
            if ~exist(fin_fmri,'file')
                preproc_fmri(mni_fmri,{'nuis'},'',1,0,motpar_fn);
            end
        end
    end
end

