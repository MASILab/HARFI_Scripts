function [] = preproc_alignment_SG(epi_dset_fn,T1_reg_dir,DISP)
% calculates alignment EPI -> T1 , T1 -> MNI, and EPI -> MNI
%
% inputs:
% ------------
% * epi_dset_fn: filename of epi dataset to register
% * T1_reg_dir: folder resulting from preproc_anat.m. Contains
% T1_bet, T1_orig, etc. Usually at <subj>/func/T1_reg/
%
% sarah 3.11.2020 based off catie's script
% ------------------------------------------------------ %
% setup
% ------------------------------------------------------ %
MNI_T1_template = '/home-local/INN_final_project/preproc/MNI152_T1_2mm_brain.nii.gz';
assert(exist(MNI_T1_template)~=0);

% check input files
if ~exist(T1_reg_dir,'dir')
    error('no reg dir: need to run preproc_anat() first');
end
if ~exist(epi_dset_fn,'file')
    error('input EPI dataset does not exist');
end

% get input EPI dimensions
hh = load_untouch_header_only(epi_dset_fn);
nslc = hh.dime.dim(4);
nframes = hh.dime.dim(5);
TR = hh.dime.pixdim(5);

% ------------------------------------------------------ %
% get path info
% ------------------------------------------------------ %
[fp,fn,ext] = fileparts(epi_dset_fn);
% make output reg dir
REG_dir = [fp,'/REG/']; % reg directory
if ~exist(REG_dir,'dir'); mkdir(REG_dir); end;
%dset_base = [REG_dir,fn];
% gunzip if necessary
if strcmp(ext(end-2:end),'.gz')
    display('gunzipping file...');
    gunzip(epi_dset_fn);
end

% ------------------------------------------------------ %
% extract reference EPI volume: middle volume of input volume (if
% one doesn't already exist. Writes as "oneVol" within REG_* dir 
% ------------------------------------------------------ %
%if ~exist([out_dir,'/oneVol.nii'],'file')
display('extracting one EPI reference vol ');
display('----------------------------------------------');

voln = floor((nframes-7)/2);

out_fn = [REG_dir,'/oneVol'];
cmd = ['3dcalc -a ', ...
       epi_dset_fn, '''[',num2str(voln),']''', ...
       ' -expr ','''a''', ...
       ' -prefix ',out_fn,'.nii'];
display(cmd); system(cmd);

% ------------------------------------------------------ %
% correct inhomogeneity in reference EPI vol
% ------------------------------------------------------ %
display('bias-field correcting EPI reference vol ');
display('----------------------------------------------');
cmd = ['fast --nopve --type=2 -B -b ', ...
       ' --out=',REG_dir,'/oneVol_HC.nii', ...
       ' ',REG_dir,'/oneVol.nii'];
display(cmd); system(cmd); 
% output is *_HC_restore.nii.
% also apply automask to this image:
cmd = ['3dAutomask -dilate 2 ', ...
       ' -apply_prefix ',REG_dir,'/oneVol_HC_am.nii', ...
       ' ',REG_dir,'oneVol_HC_restore.nii'];
display(cmd); system(cmd);
% show it
if (DISP)
    ntmp = load_untouch_nii([REG_dir,'oneVol_HC_am.nii']);
    fig; slmontage_nih(ntmp.img,[]);
    colormap gray;
end
    
% ------------------------------------------------------ %
%
% Registration to T1
%
% ------------------------------------------------------ %
% note, only rigid body is recommended for within-subject EPI to
% T1.

display('linking T1 & MNI files to "REG" dir ');
display('----------------------------------------------');
T1_bet = [T1_reg_dir,'/T1_bet.nii'];
T1_img = [T1_reg_dir,'/T1_orig.nii'];
assert(exist(T1_img)~=0)
assert(exist(T1_bet)~=0)

cmd = ['ln -s ',T1_bet,' ',REG_dir,'/T1_bet.nii']; system(cmd)
cmd = ['cp ',MNI_T1_template,' ',REG_dir,'/MNI152_T1_2mm_brain.nii.gz']; system(cmd)
cmd = ['3drefit -view orig -space ORIG ',REG_dir,'/MNI152_T1_2mm_brain.nii.gz'];

% ------------------------------------------------------ %
% Register ref vol from EPI -> TI. FSL's epi_reg works well!
% ------------------------------------------------------ %
display('*calculating* alignment of EPI to subject T1 ');
display('----------------------------------------------');

% epi reference
epi_ref = [REG_dir,'/oneVol_HC_am'];
    
% reg dir
%[fp,fn,ext] = fileparts(dset_fn); % original input                                
%reg_dir = [fp,'/reg/']; %<- dir where anat reg params are kept       
% this is now T1_reg_dir
    
% string for full epi_reg command
cmd_stem = ['epi_reg ', ...
            ' --epi=',epi_ref,'.nii', ...
            ' --t1=',T1_img, ...
            ' --t1brain=',T1_bet, ...
            ' --out=',epi_ref,'_EPI2T1.nii'];

% if we have a WM seg for this session, include in cmd to skip FAST
if exist([T1_reg_dir,'/T1_bet_fast_wmseg.nii.gz'])
    cmd = [cmd_stem, ...
           ' --wmseg=',T1_reg_dir,'/T1_bet_fast_wmseg.nii.gz'];
    display('Using existing wmseg...');
    display(cmd); system(cmd);
else
    % otherwise, run all of epi_reg
    display(cmd_stem); system(cmd_stem);
    % and copy images over for next time
    cmd = ['cp ', epi_ref,'_EPI2T1_fast_wmseg.nii.gz', ...
           ' ',T1_reg_dir,'/T1_bet_fast_wmseg.nii.gz'];
    system(cmd);
    cmd = ['cp ', epi_ref,'_EPI2T1_fast_wmedge.nii.gz', ...
           ' ',T1_reg_dir,'/T1_bet_fast_wmedge.nii.gz'];
    system(cmd);
end

% ------------------------------------------------------ %
% Apply T1 -> MNI transform to the refvol
% ------------------------------------------------------ %
% T1 to MNI
    in_fn = T1_bet;  % skull-stripped (MNI is skull-stripped)
    out_fn = [REG_dir,'reg_toMNI_T1'];
    cmd = ['flirt -in ' in_fn, ...
           ' -ref ' MNI_T1_template, ...
           ' -omat ' REG_dir,'reg_subT1_to_maMNI.mat' ...
           ' -dof 12 -out ',out_fn,'.nii'];
    display(cmd); system(cmd);
% ------------------------------------------------------ %
% Apply EPI -> MNI transform
% ------------------------------------------------------ %
% EPI to MNI  
    in_fn = [epi_ref, '_EPI2T1.nii'];  % original input
    out_fn = [REG_dir,'reg_toMNI_EPI'];
    cmd = ['flirt -in ' in_fn, ...
           ' -ref ' MNI_T1_template, ...
           ' -applyxfm -init ' REG_dir,'reg_subT1_to_maMNI.mat' ...
           ' -interp trilinear -out ',out_fn,'.nii'];
    display(cmd); system(cmd);

% ------------------------------------------------------ %
% clean up intermediate files
% ------------------------------------------------------ %
mkdir([REG_dir,'/intermed_files']);
cmd = ['mv', ...
       ' ',REG_dir,'oneVol_HC_bias.nii.gz', ...
       ' ',REG_dir,'oneVol_HC_seg.nii.gz', ...
       ' ',REG_dir,'oneVol_HC_restore.nii.gz', ...
       ' ',REG_dir,'oneVol_HC_am_EPI2T1_fast_wmseg.nii.gz', ...
       ' ',REG_dir,'oneVol_HC_am_EPI2T1_fast_wmedge.nii.gz', ...
       ' ',REG_dir,'/intermed_files/'];
display(cmd); system(cmd);

