% -------------------------------------------------------
%
%    pp12Main - loading nifti data and executing first preprocessing steps
%
%    Ver. 1.0.0
%    Created:         Daniela Zoeller      (30.10.2015)
%    Based on: preprocROItimeCourseExample.m (Jonas Richiardi and Giulia Preti)
%    Last modified:   Celine Provins      (03.2020)

%
% ------------------------------------------------------
%
% should be executed with matlab2015 and spm12
%
% images should already be converted to NIFTI
%
% 
clear all; close all;
warning('off')

%% set up path
RestingState_mode = 1; %1 = resting state data, 0 = working memory data
% celinepath = 'J:\Anjali_Diffusion_Pipeline\Celine';
celinepath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
if RestingState_mode
    datapath = fullfile(celinepath,'data','RestingState');
else
    datapath = fullfile(celinepath,'data','ControlsWM');
end
codeBasePath=fullfile(celinepath,'Preprocessing-Hub');
spmPath='/media/miplab-nas2/Data/Anjali_UB50_Data/fMRI-Sleep-Whole-Scan/Preprocess/spm12';
AALfile=fullfile(celinepath,'Preprocessing-Hub','AAL90_correctLR.nii');

addpath(genpath(fullfile(codeBasePath,'functions')));
addpath(genpath(spmPath));

%% Finding each patient folder
dirs = dir(fullfile(datapath, 's*'));
folders = cell(size(dirs));
for i = 1:size(dirs,1)
    folders{i} = fullfile(datapath,dirs(i).name);
end

%% convert 4D nifti volumes to 3D
%only for functional volumes
for j = 1:size(folders)
    cd(fullfile(folders{j},'func'));
    files = dir('*.nii');
    for i = 1 : size(files)
        spm('defaults','fmri');
        spm_jobman('initcfg');
        matlabbatch{1}.spm.util.split.vol = {fullfile(folders{j},'func',files(i).name)}; 
        matlabbatch{1}.spm.util.split.outdir = {fullfile(folders{j},'func')};
        spm_jobman('run',matlabbatch);
        delete(files(i).name);
    end
    mkdir(fullfile(folders{j},'func'),'realigned');
end
%% actual preprocessing

for j = 1:length(folders)
    fprintf('%s%d%s%d\n','Processing subject number ', j, '. Total number of subjects is ', size(folders,1));
    
    structPath = fullfile(folders{j},'struct');
    
    functPath = fullfile(folders{j}, 'func');

    if(exist(structPath,'dir') && exist(functPath,'dir'))
        if RestingState_mode
            preprocess12(functPath,structPath,'procChain',{'QC','coregReslice'},...
                    'QCcoef',1.5,'coregDirection','FtoS','atlasFile',AALfile,'atlasType','AAL');
        else
            preprocess12(functPath,structPath,'procChain',{'reset','realign','QC','coregReslice'},...
                    'QCcoef',1.5,'coregDirection','FtoS','atlasFile',AALfile,'atlasType','AAL');
        end
    
    end
end
cd(codeBasePath);

