%% Convert segmented masks to MNI space
clear all; close all;
warning('off');

%% Set up paths
dataset = 'WorkingMemory';
%celinePath = 'J:\Anjali_Diffusion_Pipeline\Celine';
celinePath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
codeBasePath = fullfile(celinePath,'kmeans');
spmPath = fullfile(celinePath,'spm12');
dataPath = fullfile(celinePath,'data',dataset);
capPath = fullfile(celinePath,'data','CAP_RS');
savePath = fullfile(celinePath,'analysis_after_clustering','masks');

addpath(genpath(spmPath));
addpath(genpath(codeBasePath));

%% Selecting patient folders
dirs = dir(fullfile(dataPath, 's*'));
subjects = cell(size(dirs));
for i = 1:size(dirs,1)
    subjects{i} = dirs(i).name;
end

%% Convert masks
cap_model = fullfile(capPath,'GSR_CAP_1_mean.nii');
cap_info = spm_vol(cap_model);
dim = cap_info.dim;
GM_mask=zeros(length(subjects),dim(1)*dim(2)*dim(3));
WM_mask=zeros(length(subjects),dim(1)*dim(2)*dim(3));
directory = pwd;

for i = 1:length(subjects)
    fprintf('%s%d%s%d\n','Processing subject number ', i, '. Total number of subjects is ', size(subjects,1));
    subject = subjects{i};
    segPath = fullfile(dataPath,subject,'struct','Segmented');
    %Normalize the masks
    cd(segPath);
    deformation = dir('y*');
    defld = fullfile(deformation.folder,deformation.name);
    NormalizeToMD(segPath,'c2',defld);
    NormalizeToMD(segPath,'c1',defld);
    
%     %redimension the mask to match dimensions of CAPs
%     c2=dir('MNI_c2*');
%     c1=dir('MNI_c1*');
%     c2_mask = mapVolumeToVolume(fullfile(c2.folder,c2.name),cap_model);
%     c1_mask = mapVolumeToVolume(fullfile(c1.folder,c1.name),cap_model);
%     %reshape to 1D
%     WM_mask(i,:) = c2_mask(:);
%     GM_mask(i,:) = c1_mask(:);
end
% keep in the mask only pixels present in 70% of the frames
% thresh = 0.7;
% WM_mask = column_nul(WM_mask,thresh);
% GM_mask = column_nul(GM_mask,thresh);
% niftiwrite(WM_mask,fullfile(savePath,'testWM.nii'));
% save(fullfile(savePath,strcat(dataset,'_masks.mat')),'WM_mask','GM_mask');
% cd(directory)
% 
% function out = column_nul(mask,thresh)
%     idx_zero = find(mean(mask,2)<=thresh);
%     for j = 1:length(idx_zero)
%         mask(:,idx_zero(j))=0;
%     end
%     out = mask;
% end
% 

   
    
 
    