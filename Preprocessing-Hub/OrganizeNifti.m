% Reorganizing the emplacement of the nifti files
% Created by Celine Provins on 25.02.2020 (celine.provins@epfl.ch)

clear all
close all
clc

%% set up path
inputPath = '/media/miplab-nas2/Data/AgCC/WorkingMemory/Controls';
outputPath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine/data/ControlsWM';
% inputPath = 'J:\AgCC\WorkingMemory\Controls';
% outputPath = 'J:\Anjali_Diffusion_Pipeline\Celine\data\ControlsWM';

%% Move all nifti files to a new folder
% Working_Memory .nii file to 'func' folder
% T1 .nii file to 'struct' folder

%Finding each patient folder
dirs = dir(fullfile(inputPath, '*'));
folders = cell(size(dirs));
for i = 1:size(dirs,1)
    folders{i} = strcat('s',dirs(i).name);
end

for j = 4:size(folders)
    fprintf('%s%d%s%d\n','Moving folders of subject number ', j, '. Total number of subjects is ', size(folders,1));
    if ~exist(fullfile(outputPath, folders{j}),'dir')
      mkdir(outputPath,folders{j});
    end
    subjectPath = fullfile(outputPath,folders{j});
    
    %moving the files
    if ~exist(fullfile(subjectPath,'func'), 'dir')
      mkdir(subjectPath,'func');
    end
    if ~exist(fullfile(subjectPath,'struct'), 'dir')
        mkdir(subjectPath,'struct');
        mkdir(fullfile(subjectPath,'struct'),'Segmented');
    end
    funcPath = fullfile(subjectPath,'func');
    structPath = fullfile(subjectPath,'struct');
    segPath = fullfile(subjectPath,'struct','Segmented');
    
%     cd(fullfile(inputPath,folders{j}(2:end),'T1'));
%     files = dir('*.nii*');
%     for i=1:size(files)
%         movefile(fullfile(inputPath,folders{j}(2:end),'T1',files(i).name),...
%             fullfile(outputPath,folders{j},'struct'));
%     end
    cd(fullfile(inputPath,folders{j}(2:end),'Working_Memory'));
    files = dir('*.nii*');
    for i = 1:length(files)
            copyfile(fullfile(files(i).folder,files(i).name),funcPath);
    end
    
    cd(fullfile('/media/miplab-nas2/Data/AgCC/Anjali/Controls',folders{j},'struct'))
    rs = dir('rs*.nii');
    c1 = dir('c1rs*.nii');
    c2 = dir('c2rs*.nii');
    c3 = dir('c3rs*.nii');
    iy = dir('iy_rs*.nii');
    y = dir('y_rs*.nii');
    copier(rs,structPath);
    copier(c1,segPath);
    copier(c2,segPath);
    copier(c3,segPath);
    copier(iy,segPath);
    copier(y,segPath);
    
    %gunzip .nii.gz files to extract nifti
    fprintf('%s%d%s%d\n','Gunzip for subject number ', j, '. Total number of subjects is ', size(folders,1));
    cd(fullfile(outputPath,folders{j},'func'));
    files = dir('*.gz');
    for i = 1:size(files)
        gunzip(files(i).name);
        delete(files(i).name);
    end
end
cd('/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine/Preprocessing-Hub');

function copier(file,outputPath)
    copyfile(fullfile(file.folder,file.name),outputPath);
end