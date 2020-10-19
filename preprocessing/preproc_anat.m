function [] = preproc_anat_SG(anat_fn, MNI_T1_template)
% function [] = preproc_anat(anat_fn)
% pre-process T1 image: performs skull-strip (via BET) and

% version 11 March 2020
% sarah based on catie's script
%
% INPUT: anat_fn: full path to T1 MPRAGE (nii format)
%% defaults:
% MNI_T1_template = '/home-local/INN_final_project/preproc/MNI152_T1_2mm_brain.nii.gz';


[fp,fn,ext] = fileparts(anat_fn);
out_dir = [fp,'/T1_reg/']; 
display('----------------------------------------------');
display(['output of anatomic bet and registration will be in: ']);
display(out_dir);
display('----------------------------------------------');
if ~exist(out_dir)
    mkdir(out_dir);
end

% ------------------------------------------------------ %
% add symlink to original T1 image
% ------------------------------------------------------ %
cmd = ['ln -s ', ...
       anat_fn, ...
       ' ',out_dir,'/T1_orig.nii'];
display(cmd); system(cmd);
% ------------------------------------------------------ %
%
% add symlink to MNI template (actually, copy it)
% ------------------------------------------------------ %
cmd = ['cp ', ...
       MNI_T1_template, ...
       ' ',out_dir];
system(cmd);
% change it to a "+orig" 
[~,MNInm,MNIex] = fileparts(MNI_T1_template);
cmd = ['3drefit -view orig -space ORIG ',out_dir,MNInm,MNIex];
display(cmd); system(cmd);

% ------------------------------------------------------ %
%
% try skull-strip of our T1 here
%  preproc_fmri.m assumes output is reg/T1_bet.nii
% ------------------------------------------------------ %
% this is highly recommended, both for t1->mni and epi->t1
display('running BET on anatomic... ');
display('----------------------------------------------');
cmd = ['bet ',anat_fn, ...
       ' ',out_dir,'T1_bet.nii', ...
       ' -R -m -o -f 0.4'];
display(cmd); system(cmd);
cmd = ['gunzip ',out_dir,'/T1_bet.nii.gz'];
system(cmd);

