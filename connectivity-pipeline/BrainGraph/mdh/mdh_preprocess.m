




function  G = mdh_preprocess(subjID,hcp_root) % Based on sgwt_preprocess.m
% MDH_PREPROCESS Creates necessary files for building GM, WM,GM+WM graphs.
%
%
%
% Contributers: 
%    Martin Larsson
%    David Abramian
%    Hamid Behjat
% 2017-2018.


% Code revised by: Anjali Tarun:
% Changes made:
%     -- to use whole brain volume instead of WM only
%     -- obtains whole brain brain mask
% February 2018

f = struct;

% Folders.
f.T1w = fullfile(hcp_root,subjID,'T1w');
f.diffusion = fullfile(f.T1w,'Diffusion');
f.extracted = fullfile(f.T1w,'extracted');
f.graph = fullfile(f.T1w,'graph');
if ~exist(f.extracted,'dir')
    mkdir(f.extracted);
end
if ~exist(f.graph,'dir')
    mkdir(f.graph);
end

% Files.
f.ref_125 = fullfile(f.diffusion,'nodif_brain_mask.nii');
f.aparc_aseg_07 = fullfile(f.T1w,'aparc+aseg.nii');
f.aparc_aseg_125 = fullfile(f.T1w,'aparc+aseg_1.25.nii');
f.ribbon_125 = fullfile(f.extracted,'ribbon_1.25.nii');
f.white_125 = fullfile(f.extracted,'white_1.25.nii');
f.brainmask_fs = fullfile(f.T1w, 'brainmask_fs.nii');
f.brainmask_fs_125 = fullfile(f.T1w, 'brainmask_fs_125.nii');
f.corpus_callosum_125 = fullfile(f.extracted,'corpus_callosum_1.25.nii');
f.ribbon_125_conn = fullfile(f.extracted,'ribbon_1.25_connected.nii');
f.gmwm_125 = fullfile(f.extracted,'gmwm_125.nii');  
f.T1_125 = fullfile(f.T1w, 'T1w_acpc_dc_restore_1.25.nii');

% extracts 
% Downsample aparc+aseg.
ml_extract(f.aparc_aseg_07);
ml_extract(f.ref_125);
if ~exist(f.aparc_aseg_125,'file')
    f.aparc_aseg_125 = da_resize(f.aparc_aseg_07,f.ref_125);
end

if ~exist(f.aparc_aseg_125,'file')
    f.aparc_aseg_125 = da_resize(f.aparc_aseg_07,f.ref_125);
end

if ~exist(f.brainmask_fs_125,'file')
    ml_extract(f.brainmask_fs);
    
    f.brainmask_fs_125 = da_resize(f.brainmask_fs,f.ref_125);
    % Clean mask.
    [wh,wv] = ml_load_nifti(f.brainmask_fs_125);
    cc = bwconncomp(wv,26);
    [~,ind] = max(cellfun(@length,cc.PixelIdxList));
    wv = zeros(size(wv));
    wv(cc.PixelIdxList{ind}) = 1;
    spm_write_vol(wh,wv);
end


% Extract volumes.
ids_ribbon = [3 42 1000:2999];
ids_corpus_callosum = [86 251:255];
ids_white = [2 41 ids_corpus_callosum];

if ~exist(f.corpus_callosum_125,'file')
    ml_extract_volumes(f.aparc_aseg_125,f.corpus_callosum_125,...
        ids_corpus_callosum);
end

if ~exist(f.ribbon_125,'file')
    ml_extract_volumes(f.aparc_aseg_125,f.ribbon_125,ids_ribbon);
end

if ~exist(f.white_125,'file')
    ml_extract_volumes(f.aparc_aseg_125,f.white_125,ids_white);
    
    % Clean white matter.
    [wh,wv] = ml_load_nifti(f.white_125);
    cc = bwconncomp(wv,26);
    [~,ind] = max(cellfun(@length,cc.PixelIdxList));
    wv = zeros(size(wv));
    wv(cc.PixelIdxList{ind}) = 1;
    spm_write_vol(wh,wv);
end



G = struct;
G.f = f;
