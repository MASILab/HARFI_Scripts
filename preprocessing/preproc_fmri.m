function OUT = preproc_fmri_SG(varargin);
%function [] = preproc_fmri(dset_fn, CHAIN, motpars, DISP)
% pre-process fmri data
% *  first we run dicoms2nii.py to create nii for each
% *  CHAIN is the order of preprocessing steps to do. Cell array:
% {'drop','ret','motco','sltim','nuis','ssm','reg_setup','align_T1','align_MNI'};
% note.. will not change the order from above, even if you pass in
% a re-ordered array.
% *  if we start from a later stage, DON'T pass later input to
% dset_fn. Will check for existing files.
% * pass in motion-coreg alignment mats (*.volreg_mats.aff12.1D) if we want to use existing motion pars for realignment.
% * DISP to display QA plots
%
% @author catie  8/8/17 
% version: 9/24/17, moved the T1 alignment parts to
% preproc_alignment.m. Now, CHAIN may contain only:  {'drop','ret','motco','sltim','nuis','ssm'} 

OUT=[];
% ------------------------------------------------------ %
% args
% ------------------------------------------------------ %
if length(varargin)<2
    error('need to specify dset & CHAIN');
end
dset_fn = varargin{1};
CHAIN = varargin{2};
volreg_mats = [];
DISP = 1;
if length(varargin)>=3
    volreg_mats = varargin{3};
end
if length(varargin)==4
    DISP = varargin{4};
end

% ------------------------------------------------------ %
% get path info
% ------------------------------------------------------ %
[fp,fn,ext] = fileparts(dset_fn);
% make output dir
% each scan will have its own dir in "func/proc"
out_dir = [fp,'/proc/',fn,'/']; 
if ~exist(out_dir,'dir'); mkdir(out_dir); end;
dset_base = [out_dir,fn];

display(['*********************************************']);
display(['input file: ',dset_fn]);
display(['*********************************************']);
% ------------------------------------------------------ %
% get dataset info
% ------------------------------------------------------ %
hh = load_untouch_header_only(dset_fn);
nslc = hh.dime.dim(4);
nframes = hh.dime.dim(5);
TR = hh.dime.pixdim(5);

display(['read dimensions: nframes=',num2str(nframes), ...
         ' , nslc=',num2str(nslc), ', TR=',num2str(TR)])
display('assuming Siemens interleaved up slice order....');
if (mod(nslc,2)==0) 
    % Siemens convention for interleaved even #'d slices
    slice_order = [2:2:nslc,1:2:nslc-1];
    tpattern = 'alt+z2';
else
    % Siemens convention for interleaved odd #'d slices
    slice_order = [1:2:nslc,1:2:nslc-1];
    tpattern = 'alt+z';
end
display(num2str(slice_order));


% ------------------------------------------------------ %
% initialize the current output fn
% ------------------------------------------------------ %
curr_out_fn = [dset_base]; % this is updated each time we run a new
                           % preproc step

% ------------------------------------------------------ %
%
% drop initial time frames 
%
% ------------------------------------------------------ %
if (ismember('drop',CHAIN))
    display('dropping initial frames ');
    display('----------------------------------------------');
    in_fn = dset_fn;  % original input
    out_fn = [dset_base,'_dr'];
    cmd = ['3dcalc -a ', ...
           in_fn, '''[7..$]''', ...
           ' -expr ','''a''', ...
           ' -prefix ',out_fn,'.nii'];
    display(cmd); system(cmd);
    curr_out_fn = out_fn; % reset current out
end
% ------------------------------------------------------ %
%
% RETROICOR
%
% ------------------------------------------------------ %
if (ismember('ret',CHAIN))
    display('retroicor ');
    display('----------------------------------------------');
    if ~exist(phys_fn,'file')
        display('*** WARNING, no phys file - continuing ***');
    else
        physio_mat = load(phys_fn); %<- load an output name from physio correction.
    
        in_fn = curr_out_fn;             %[dset_base,'_dr'];
        out_fn = [curr_out_fn,'_re'];    %[dset_base,'_dr_re'];
        nn = load_untouch_nii([in_fn,'.nii']);
        image_matrix = nn.img;
        
        % run retroicor without RVHR
        [image_matrix_corrected,PHASES,REGRESSORS,OTHER] = ...
            retroicor_rvhr(image_matrix,slice_order,TR, ...
                           physio_mat.card_times, ...
                           physio_mat.resp,[]);
        
        nnw = nn;
        nnw.img = image_matrix_corrected;
        save_untouch_nii(nnw,[out_fn,'.nii']);
        curr_out_fn = out_fn; % reset current out
    end
end
% ------------------------------------------------------ %
%
% motion coregistration
%
% ------------------------------------------------------ %
% alas, motion params *.mat/MAT_ are written with one fewer zero
% than applyxfm4D expects. for now, use 3dvolreg..

if (ismember('motco',CHAIN))
    display('motion coreg ');
    display('----------------------------------------------');
  
    in_fn = curr_out_fn          %[dset_base,'_dr_re'];
    out_fn = [curr_out_fn,'_mo'] % [('volreg') dset_base,'_dr_re_mc'];
    
    if isempty(volreg_mats)  % run coreg from scratch
                             % ---------------------------
        
        % cmd = ['mcflirt -in ',in_fn, ...
        %       ' -out ',out_fn,'.nii', ...
        %       ' -mats ', ...
        %       ' -stages 3 -plots'];
        %display(cmd);  system(cmd);
        %curr_out_fn = out_fn; % reset current out
        %motpar_fn = [out_fn,'.nii.par']; 
        
        % the *aff12.1D is what we apply to new one.
        cmd = ['3dvolreg -verbose ', ...
              ' -1Dmatrix_save ',in_fn,'.volreg_mats', ...
              ' -1Dfile ',in_fn,'.volreg_par', ...
              ' -base ', in_fn,'.nii''[0]''', ...
              ' -prefix ',out_fn,'.nii', ...
              ' ',in_fn,'.nii'];
        system(cmd);
        motpar_fn = [in_fn,'.volreg_par'];
	motpar_mtx = [in_fn,'.volreg_mats.aff12.1D'];
    else  % if we want to just apply the existing pars:
          % ----------------------------------------------        
        display(['motion coreg using existing file:  ',volreg_mats]);
        % should be: *'.volreg_mat1D.aff12.1D'
        cmd = ['3dAllineate ', ...
               ' -1Dmatrix_apply ',volreg_mats, ... 
               ' -base ',in_fn,'.nii''[0]''', ...
               ' -master ',in_fn,'.nii', ...
               ' -prefix ',out_fn,'.nii', ...
               ' -source ',in_fn,'.nii'];
        display(cmd); 
        system(cmd);
                     
        tmp = findstr(volreg_mats,'_mats');            
        motpar_fn = [volreg_mats(1:tmp-1),'_par'];
        motpar_mtx = volreg_mats;
	% yes, works!
        
        %cmd = ['applyxfm4D ', ...
        %       ' ',in_fn,'.nii', ... 
        %       ' ',ref_vol, ...
        %       ' ',out_fn, ...
        %       ' ',motcor_mats];
        %motpar_fn = [motcor_mats(1:end-4),'.par']; % for later QA.

    end
    curr_out_fn = out_fn; % reset current out
    OUT.motpar_mtx = motpar_mtx; % affine xform          
    OUT.motpar_fn = motpar_fn; % 6 dof ts (for regression)
end
    
% ------------------------------------------------------ %
%
% slice-timing correction
%
% ------------------------------------------------------ %
if (ismember('sltim',CHAIN))
    display('slice-timing correction ');
    display('----------------------------------------------');
    in_fn = curr_out_fn;         %[dset_base,'_dr_re_mc'];
    out_fn = [curr_out_fn,'_st']; %[dset_base,'_dr_re_mc_st'];
    cmd = ['3dTshift ', ...
           ' -tpattern ',tpattern, ...   
           ' -Fourier ', ...
           ' -prefix ',out_fn,'.nii', ...
           ' ',in_fn,'.nii'];
    display(cmd);  system(cmd);
    curr_out_fn = out_fn; % reset current out
end
% ------------------------------------------------------ %
%
% nuisance regression
% (later do with & without RV/HR)
%
% ------------------------------------------------------ %
if (ismember('nuis',CHAIN))
    display('regressing motion parameters & trends ');
    display('----------------------------------------------');
    % only regress motion parameters & low-order polynomials
    % (4th-order for 15min scans)
    in_fn = curr_out_fn;              %[dset_base,'_dr_re_mc_st'];
    out_fn_tmp = [curr_out_fn,'_nr_tmp']; %[dset_base,'_dr_re_mc_st_nr_tmp'];
    out_fn = [curr_out_fn,'_nr'];
    % check for motion parameters
    if ~exist(motpar_fn,'file')
        error('no motion params!');
    end
    % detrend 
    cmd = ['3dDetrend ', ...
           ' -prefix ',out_fn_tmp,'.nii', ...
           ' -polort 4', ...
           ' -vector ',motpar_fn, ...
           ' ',in_fn,'.nii'];
    display(cmd); system(cmd); 
    %% add back the mean of original 
    % calculate mean of the input image
    cmd = ['3dTstat -mean ', ...
           ' -prefix ',out_dir,'/meanvol_tmp.nii', ...
           ' ',in_fn,'.nii'];
    display(cmd); system(cmd);
    cmd = ['3dcalc -a ', out_fn_tmp,'.nii', ...
           ' -b ', out_dir,'/meanvol_tmp.nii', ...
           ' -expr ''a+b''', ...
           ' -prefix ', out_fn,'.nii'];
    display(cmd); system(cmd);
    % reset current out
    curr_out_fn = out_fn;
end


% ------------------------------------------------------ %
% extract reference EPI volume: middle volume of mot-corr data
% and automask it
% ------------------------------------------------------ %
% do this if a motion-coreg was specified in our chain
% and we don't already have one
if (ismember('motco',CHAIN) && ~exist([out_dir,'/oneVol.nii'],'file'))
    display('extracting one EPI reference vol ');
    display('----------------------------------------------');

    voln = floor((nframes-7)/2);

    % locate & gunzip the *mo file
    motcor_fn = [];
    dd = dir(out_dir); 
    for kk=1:length(dd)
        if length(strfind(dd(kk).name,'_mo.nii'))>0
            motcor_fn = [dd(kk).folder,filesep,dd(kk).name];
            break;
        end
    end
    if isempty(motcor_fn)
        error('no motion-coreg dataset found');
    elseif strcmp(motcor_fn(end-2:end),'.gz')
        display('gunzipping mot-corr dataset...');
        gunzip(motcor_fn);
        motcor_fn = motcor_fn(1:end-3);
    end
    % extract the single volume
    out_fn = [out_dir,'/oneVol'];
    cmd = ['3dcalc -a ', ...
           motcor_fn, '''[',num2str(voln),']''', ...
           ' -expr ','''a''', ...
           ' -prefix ',out_fn,'.nii'];
    display(cmd); system(cmd);
end
% automask 
% ------------------------------------------------------ %
if ~exist([out_dir,'/autoMask.nii'],'file');
    display('running automask ');
    display('----------------------------------------------');
    in_fn = [out_dir,'/oneVol'];  % image on which to calculate mask
    cmd = ['3dAutomask ', ...
           ' -dilate 2', ...
           ' -prefix ',out_dir,'/autoMask.nii', ...
           ' ',in_fn,'.nii'];
    display(cmd); system(cmd);
    
    % apply this mask to our (latest) dataset:
    %cmd = ['3dcalc -prefix ',curr_out_fn,'.nii', ...
    %       ' -a ',dset_base,'_dr_re_mc_st_nr.nii', ...
    %       ' -b ',out_dir,'/autoMask.nii', ...
    %       ' -expr ''a*b'' '];
    %display(cmd); system(cmd);
end

% ------------------------------------------------------ %
%
% spatial blurring (native space)
%
% ------------------------------------------------------ %
if (ismember('ssm',CHAIN))
    display('running spatial smooth (4mm) ');
    display('----------------------------------------------');
    FWHM = 4; % mm
    
    in_fn = curr_out_fn;               %[dset_base,'_dr_re_mc_st_nr'];
    out_fn = [curr_out_fn,'_sm'];      %[dset_base,'_dr_re_mc_st_nr_sm'];
    if ~exist([out_fn,'.nii'],'file')
        cmd = ['3dBlurToFWHM ', ...
               ' -mask ',out_dir,'/autoMask.nii', ...
               ' -input ' in_fn,'.nii', ...
               ' -prefix ' out_fn,'.nii', ...
               ' -FWHM ',num2str(FWHM)];
        display(cmd); system(cmd);
    else
       display('...already done.'); 
    end
end % don't update curr_out_fn this time.
    

% ---------------------------------
%
% some QA plots
% global signal, tsnr and mean and sd maps 
%
% --------------------------------
if (ismember('qa',CHAIN))
    QA_dir = [out_dir,'/QA_files/'];
    mkdir(QA_dir);

    % ---------------------------------------------------------------- %
    % motion params
    if exist(motpar_fn,'file')
        motpar = load(motpar_fn);
        figure(5); set(gcf,'color','w');
        plot(motpar); ylabel('mm'); xlabel('frame');
        set(gcf,'position',[14   237   770   333]);
        set(gcf,'paperpositionmode','auto');
        title(dset_fn)
        print(gcf,'-depsc',[QA_dir,'/QA_motpars']);
    end
    % ---------------------------------------------------------------- %


    % ----------------------------------------------------------------
    % 

    % global signal & GS corr  (should be detrended)
    in_fn = curr_out_fn ; % latest 
    nn = load_untouch_nii([in_fn,'.nii']);
    mm = load_untouch_nii([out_dir,'/autoMask.nii']);
    tmp = size(nn.img); dimvec = tmp(1:3);
    assert(all(dimvec==size(mm.img)))

    mask = find(mm.img>0);
    V = vol_4d2d(nn.img); 
    Y = V(mask,:)';
    gs = mean(Y,2); 
    gscorr = corr(gs,double(Y));
    gsvol = zeros(dimvec); gsvol(mask)=gscorr;

    % save plots
    figure(2); set(gcf,'color','w');
    rowPlots_nih({gsvol},4:2:27, [], {}, '');
    set(gcf,'position',[31   123   932   108]);
    set(gcf,'paperpositionmode','auto');
    print(gcf,'-depsc',[QA_dir,'/QA_gscorr']);

    figure(3); set(gcf,'color','w');
    plot(gs); title('global signal');
    set(gcf,'position',[40   106   939   167]);
    set(gcf,'paperpositionmode','auto');
    print(gcf,'-depsc',[QA_dir,'/QA_gssig']);
    % ---------------------------------------------------------------- %


    % ---------------------------------------------------------------- %
    % compute stdev & tsnr map on latest image (should be detrended)

    % find the "nv" filename:
    nr_fn = [];
    dd = dir(out_dir); 
    for kk=1:length(dd)
        if length(strfind(dd(kk).name,'_nr'))>0
            nr_fn = [dd(kk).folder,filesep,dd(kk).name];
            break;
        end
    end


    if ~isempty(nr_fn)
        
        % mean map
        display('Computing mean map (*_nr.nii) .........');
        cmd = ['3dTstat	', ...
               ' -mean', ...
               ' -prefix ',QA_dir,'/meanvol.nii', ...
               ' ',nr_fn];
        display(cmd); system(cmd);
        % stdev map
        display('Computing stdev map (*_nr.nii) .........');
        cmd = ['3dTstat	', ...
               ' -stdev', ...
               ' -prefix ',QA_dir,'/stdvol.nii', ...
               ' ',nr_fn];
        display(cmd); system(cmd);
        % tSNR map
        display('Computing tSNR map (*_nr.nii) .........');
        cmd = ['3dTstat ', ...
               ' -cvarinvNOD', ...
               ' -mask ',out_dir,'/autoMask.nii', ...
               ' -prefix ',QA_dir,'/tSNR.nii', ...
               ' ',nr_fn];
        display(cmd); system(cmd);

        % display in rows: mean, sd, tsnr
        nn1 = load_untouch_nii([QA_dir,'/meanvol.nii']);
        nn2 = load_untouch_nii([QA_dir,'/stdvol.nii']);
        nn3 = load_untouch_nii([QA_dir,'/tSNR.nii']);
        clims_set = {[0 900],[0,15],[0 200]};
        figure(1); 
        rowPlots_nih({nn1.img,nn2.img,nn3.img},4:2:27, clims_set, {}, '');
        colormap gray;
        set(gcf,'position',[ 21          82        1023         366]);
        set(gcf,'paperpositionmode','auto');
        saveas(gcf,[QA_dir,'/QA_tsnrChecks.fig']);
        print(gcf,'-depsc',[QA_dir,'/QA_tsnrChecks']);
        % ---------------------------------------------------------------- %

    end
    
end % end QA






