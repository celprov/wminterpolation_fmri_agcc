clear all; close all;

%% Set up the paths
dataPath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine/data/ControlsWM';
codeBasePath='/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine/Preprocessing-Hub';
spmPath='/media/miplab-nas2/Data/Anjali_UB50_Data/fMRI-Sleep-Whole-Scan/Preprocess/spm12';
dwiPath = '/media/miplab-nas2/Data/AgCC/Anjali/Controls';

addpath(genpath(fullfile(codeBasePath,'functions')));
addpath(genpath(spmPath));

%% Finding each patient folder
dirs = dir(fullfile(dataPath, 'regFail_*'));
folders = cell(size(dirs));
for i = 1:size(dirs,1)
    folders{i} = dirs(i).name;
end

%% Coregistration
for j = 1:length(folders)
    fprintf('%s%d%s%d\n','Processing subject number ', j, '. Total number of subjects is ', size(folders,1));
    funcPath = fullfile(dataPath,folders{j},'func','realigned');
    structPath = fullfile(dataPath,folders{j},'struct');
    means = dir(fullfile(funcPath, 'mean*'));
    tmp_others=cellstr(spm_select('List',funcPath,'r*.nii'));
    tmp_others=cellfun(@(x) fullfile(funcPath,x),tmp_others,'UniformOutput',false);
    dwi = fullfile(dwiPath,folders{j}(end-3:end),'DWI','dwi_dpb_3050_norm.nii,1');

    % co-register
    spm_jobman('initcfg');
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {dwi};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {fullfile(funcPath, means(1).name)};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = tmp_others;
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'reg';
    spm_jobman('run',matlabbatch);
    clear matlabbatch
end