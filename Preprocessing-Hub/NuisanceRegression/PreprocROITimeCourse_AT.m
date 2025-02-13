%% Perform nuisance regression

%     Created by: Anjali Tarun
%     Date created: 28 March 2017
%     Modified by: Celine Provins (04.2020)
% 
%     Based on preprocROItimeCourse.m by Giulia Preti

%     Generates the files c2 and c3 with prefix 'w', needed for nuissance
%     regression and calls 'RegressOutNuissance.m' to actually perform the
%     nuisance regression.
      
%     This scripts calls functions nii_mask2target_AT.m and
%     overlapMaskNative_AT.m and requires the white matter and csf masks
%
%     Needs to be run with SPM8

clear all; close all;

%% Set up path 
dataset = 'RestingState';

%celinePath = 'J:\Anjali_Diffusion_Pipeline\Celine';
celinePath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
dataBasePath = fullfile(celinePath,'data',dataset);
codeBasePath = fullfile(celinePath,'Preprocessing-Hub');
spmPath = fullfile(codeBasePath,'spm8');
maskPath = fullfile(codeBasePath,'NuisanceRegression');
addpath(genpath(codeBasePath))
addpath(genpath(maskPath))
addpath(genpath(spmPath))

dirs = dir(fullfile(dataBasePath, 's*')); % first letter of subject folders

for i = 1:length(dirs) % iterate for each subject
    fprintf('%s%d%s%d\n','Processing subject number ', i, '. Total number of subjects is ', size(dirs,1));
    paths.r = fullfile(dataBasePath, dirs(i).name, 'func', 'realigned');
    paths.sseg= fullfile(dataBasePath, dirs(i).name, 'struct','Segmented');
    
    % Deleting the existing files
%     cd(paths.r);
%     directory2 = dir(fullfile(paths.r, 's6*'));
%     for j = 1:length(directory2)
%         delete(directory2(j).name)
%     end
%     directory2 = dir(fullfile(paths.r, 'Cov*'));
%     for j = 1:length(directory2)
%         delete(fullfile(directory2(j).folder,directory2(j).name));
%     end
    
    % warp CSF and WM to subject space
    Masks={'WhiteMask_09_121x145x121.nii','CsfMask_07_121x145x121.nii'};
    cd(paths.sseg);
    deffield=dir(fullfile(paths.sseg,'iy*.nii')); % from "new sgment", template to native space deformation field
    for j=1:length(Masks)
        tmp=fullfile(maskPath,Masks{j});
        if ~isempty(tmp) && ~exist(fullfile(paths.sseg,['w',Masks{j}]),'file')
            fprintf('Warping mask %s to native space...\n',Masks{j});
            nii_mask2target_AT(tmp,fullfile(paths.sseg,deffield.name));

        elseif isempty(tmp)
            error('No mask of right size found. Please reslice the masks of DPARSFA to the dimensions 121x145x121');
        elseif exist(fullfile(paths.sseg,['w',Masks{j}]),'file')
            fprintf('Warped mask %s already exists\n',['w',Masks{j}]);
        end
    end

    %find overlap between mask and segmentation
    overlapMaskNative_AT(paths.sseg) % saves c2 and c3 with prefix 'w'
    
    % Calls actual nuissance regression
    RegressOutNuissance(paths,0,0)
    
end
cd(maskPath);