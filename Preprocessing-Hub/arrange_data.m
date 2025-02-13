% Copying the original files to your own folder
% Created by Celine Provins 25.02.2020

clear all; close all;

inputPath = 'J:\Anjali_Diffusion_Pipeline\Celine\data\WorkingMemory';
outputPath = 'J:\Anjali_Diffusion_Pipeline\Celine\data\WorkingMemory';

%% Finding each patient folder
cd(inputPath);
dirs = dir(fullfile(inputPath, 's*'));
folders = cell(size(dirs));
for i = 1:size(dirs,1)
    folders{i} = dirs(i).name(3:end);
end

%% Copying the files from inputPath to outputPath
for j = 1:length(folders)
    if ~exist(fullfile(outputPath,folders{j}),'dir')
        % create folders in outputPath
        mkdir(outputPath,folders{j})
        subjectPath = fullfile(outputPath,folders{j});
        mkdir(subjectPath,'func');
        mkdir(fullfile(subjectPath,'func'),'realigned');
        mkdir(subjectPath,'struct');
        funcPath = fullfile(subjectPath,'func','realigned');
        structPath = fullfile(subjectPath,'struct');
        mkdir(structPath,'Segmented');
        segPath = fullfile(structPath,'Segmented');
        
        %fetching fils in inputPath and copying them in outputPath
        cd(fullfile(inputPath,['X_' folders{j}],'struct'));
        rs = dir('rs*.nii');
        cd(fullfile(inputPath,['X_' folders{j}],'struct','Segmented'));
        c1 = dir('c1rs*.nii');
        c2 = dir('c2rs*.nii');
        c3 = dir('c3rs*.nii');
        iy = dir('iy_rs*.nii');
        y = dir('y_rs*.nii');
        if ~isempty(rs)
            copier(rs,structPath);
            copier(c1,segPath);
            copier(c2,segPath);
            copier(c3,segPath);
            copier(iy,segPath);
            copier(y,segPath);
        end
        w = dir('w*.nii');
        for i = 1:length(w)
            copyfile(fullfile(w(i).folder,w(i).name),segPath);
        end
        
        cd(fullfile(inputPath,folders{j},'func','realigned'));
        mean = dir('mean*');
        rp = dir('rp*.txt');
        artifacts = dir('*artifacts.mat');
        copier(mean,funcPath);
        copier(artifacts, funcPath);
        if ~isempty(rp)
            copier(rp,funcPath);
        end
        r = dir('r*.nii');
        for i = 1:length(r)
            copyfile(fullfile(r(i).folder,r(i).name),funcPath);
        end
        cd(fullfile(inputPath,folders{j},'func'));
        nii = dir('*.nii');
        for i = 1:length(nii)
            copyfile(fullfile(nii(i).folder,nii(i).name),fullfile(subjectPath,'func'));
        end   
    end
end

function copier(file,outputPath)
    copyfile(fullfile(file.folder,file.name),outputPath);
end