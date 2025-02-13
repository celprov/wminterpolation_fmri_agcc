function nii_mask2target_AT(Mask,iyfield,intp)
% warps CSF and WM masks to native space, or viceversa to MNI (depending on
% the use of iy or y deffield)
% modelled on Chris Rorden's nii_template2target
% http://www.mccauslandcenter.sc.edu/CRNL/sw/spm8/nii_template2target.m.txt

if nargin <1 %no mask
 Mask = spm_select(1,'image','Select mni image');
end;
if nargin <2 %no field
 iyfield = spm_select(1,'image','Select ''iy'' fieldmap');
end;
if nargin <3
    intp=1;
end

[pth,nam,ext] = spm_fileparts(deblank(Mask(1,:)));
Mask = [pth,filesep,nam,ext];

[pth,nam,ext] = spm_fileparts(deblank(iyfield(1,:)));
iy = [pth,filesep,nam,ext];

spm_jobman('initcfg');
matlabbatch{1}.spm.util.defs.comp{1}.def = {iy}; % deformation field
matlabbatch{1}.spm.util.defs.ofname = '';
matlabbatch{1}.spm.util.defs.fnames = { Mask }; % image to write
matlabbatch{1}.spm.util.defs.savedir.savepwd = 1;
matlabbatch{1}.spm.util.defs.interp = intp; % degree of B-spline (from 0 to 7)
startdir = pwd;
cd(pth); %so we can write warped templates to this folder
spm_jobman('run',matlabbatch);
cd(startdir);

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
