%function nii_mask2target_AT(Mask,iyfield,intp)
% warps CSF and WM masks to native space, or viceversa to MNI (depending on
% the use of iy or y deffield)
% modelled on Chris Rorden's nii_template2target
% http://www.mccauslandcenter.sc.edu/CRNL/sw/spm8/nii_template2target.m.txt

% if nargin <1 %no mask
%  Mask = spm_select(1,'image','Select mni image');
% end;
% if nargin <2 %no field
%  iyfield = spm_select(1,'image','Select ''iy'' fieldmap');
% end;
% if nargin <3
%     intp=1;
% end

%% Set up the paths
dataset = 'RestingState';
% celinepath = 'J:\Anjali_Diffusion_Pipeline\Celine';
celinepath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
datapath = fullfile(celinepath,'data',dataset);
codeBasePath=fullfile(celinepath,'Preprocessing-Hub');
spmPath='/media/miplab-nas2/Data/Anjali_UB50_Data/fMRI-Sleep-Whole-Scan/Preprocess/spm12';
Mask = fullfile(celinepath,'Preprocessing-Hub','AAL90_correctLR.nii');

addpath(genpath(spmPath));

[pth,nam,ext] = spm_fileparts(deblank(Mask(1,:)));
Mask = [pth,filesep,nam,ext];

intp = 1;
    
dirs = dir(fullfile(datapath, 's007*'));
folders = cell(size(dirs));
for i = 1:size(dirs,1)
    folders{i} = fullfile(datapath,dirs(i).name);
end

for j = 1:length(folders)
    fprintf('%s%s\n','Normalizing atlas for subject ', folders{j});
    segPath = fullfile(folders{j},'struct','Segmented');
    cd(segPath);
    file = dir('iy*');
    iyfield = fullfile(segPath,file.name);

    [pth,nam,ext] = spm_fileparts(deblank(iyfield(1,:)));
    iy = [pth,filesep,nam,ext];

    spm_jobman('initcfg');
    matlabbatch{1}.spm.util.defs.comp{1}.def = {iy}; % deformation field
    matlabbatch{1}.spm.util.defs.ofname = '';
    matlabbatch{1}.spm.util.defs.fnames = { Mask }; % image to write
    matlabbatch{1}.spm.util.defs.savedir.savepwd = 1;
    matlabbatch{1}.spm.util.defs.interp = intp; % degree of B-spline (from 0 to 7)
    cd(pth); %so we can write warped templates to this folder
    spm_jobman('run',matlabbatch);
    cd(codeBasePath);
end

%%
% I would like to use a voxel size that is different from the one specified
% in the deformation field: 
%   myjob.comp{1}.def = {'my_deformation_field.nii'};
%   myjob.comp{2}.idbbvox.vox = [3 3 3];
%   myjob.comp{2}.idbbvox.bb = [[-93 -129 -75]; [90 90 108]]; % bounding
%   box (which portion of MNI space is incorporated, first row: start [x y z], 2nd: end; needs to be divisable by voxel size)
%   myjob.ofname = '';
%   myjob.fnames = {'my_image_to_write'};
%   myjob.savedir.saveusr = {'my_outdir'};
