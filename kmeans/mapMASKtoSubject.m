% celinePath = 'J:\Anjali_Diffusion_Pipeline\Celine'
celinePath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
dataPath = fullfile(celinePath,'data');
codePath = fullfile(celinePath,'kmeans');
spmPath = '/media/miplab-nas2/Data/Anjali_UB50_Data/fMRI-Sleep-Whole-Scan/Preprocess/spm12';
addpath(genpath(spmPath));

subjects = {'s007'};

for i = 1:length(subjects)
    fprintf('%s%s\n','Normalizing subject ', subjects{i});
    realignPath = fullfile(dataPath,subjects{i},'func','realigned');
    segPath = fullfile(dataPath,subjects{i},'struct','Segmented');
    cd(segPath);
    deformationField = dir('iy*.nii');
    cd(codePath);
    NormalizeToMD(realignPath,'s6Cov',fullfile(segPath,deformationField.name));
end