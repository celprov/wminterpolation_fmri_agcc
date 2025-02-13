%% Convert atlas to individual subject space and save them in nifti format
clear all; close all;

%% Define the paths
% celinePath = 'J:\Anjali_Diffusion_Pipeline\Celine';
celinePath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
spmPath='/media/miplab-nas2/Data/Anjali_UB50_Data/fMRI-Sleep-Whole-Scan/Preprocess/spm12';
dataPath = fullfile(celinePath,'data','ControlsRS');
Atlas = fullfile(celinePath,'analysis_after_clustering','Atlases','WM_atlas_Yeo_Combined17.nii');
addpath(genpath(spmPath));

%% Find subjects folder
dirs = dir(fullfile(dataPath, 's*'));
folders = cell(size(dirs));
for i = 1:size(dirs,1)
    folders{i} = fullfile(dataPath,dirs(i).name);
end

%% Convert the atlas for each subject
for j = 1:length(folders)
    fprintf('%s%s\n','Converting atlas for subject ', folders{j}(end-4:end));
    segPath = fullfile(folders{j},'struct','Segmented');
    
    if exist(fullfile(segPath,'wNorm_WM_atlas_Yeo_Combined17.nii'),'file')
        delete(fullfile(segPath,'wNorm_WM_atlas_Yeo_Combined17.nii'));
    end
    
    % Find all the segmented volumes of the t1 (c1, c2, c3)
    dirs = dir(fullfile(segPath, 'c*'));
    structfiles = cell(size(dirs));
    for i = 1:size(dirs,1)
        structfiles{i} = fullfile(segPath,dirs(i).name);
    end

    def = dir(fullfile(segPath,'iy*')); % inverse deformation field (iy).
    deform = fullfile(segPath, def.name);
    atlasedVol = new_Auto_Labelling(structfiles, Atlas, deform, segPath, 0.5);
end
