function WMregions = load_WM_atlas(cap_model,maskPathWM,savePath)
%--------------------------------------------------------------------------
% Created by : Celine Provins (06.2020)
%
% Load the WM atlas' regions. If exist only in nifti format, 
%   save the masks in .mat in the dimension dictated by cap_model and as a 1D vector. 
%   save the name of the WM regions
%
% INPUT
% cap_model : string containing the name of the model dictating the
%   dimension at which the masks extracted from the atlas are gonna be saved
% mastPathWM : path towards the folder containing the WM masks in nifti
%   format
% savePath : path where to save the WM regions name
%
% OUTPUT
% WMregions is a cell (1 x #regions). Each cell contains one mask as a 1D vector.
%--------------------------------------------------------------------------
cd(maskPathWM)
WMdirs = dir('Lowres*.nii');
nbr_regionWM = size(WMdirs,1);
WMregions = cell(1,nbr_regionWM);
region_nameWM = cell(1,nbr_regionWM);   
for i = 1:size(WMdirs,1)
    region_nameWM{i} = WMdirs(i).name(8:end-4);
    if i == size(WMdirs,1)
        save(fullfile(savePath,'region_nameWM.mat'),'region_nameWM');
    end
    mask_name = WMdirs(i).name(1:end-4);
    maskPath = WMdirs(i).folder;
    if ~exist(fullfile(maskPath,strcat(mask_name,'.mat')),'file')
        %Process the masks and save them in .mat format
        disp('Processing the masks')
        mask = mapVolumeToVolume(strcat(mask_name,'.nii'),cap_model);
        mask = mask(:);
        WMregions{i} = mask;
        save(fullfile(maskPath,strcat(mask_name,'.mat')),'mask');
    else
        %Load the masks
        mask_struct = load(fullfile(maskPath,strcat(mask_name,'.mat')));
        WMregions{i} = mask_struct.mask;
    end
end
