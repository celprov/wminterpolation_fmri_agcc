function regions = load_GM_atlas(savename,dirs,cap_model)
%--------------------------------------------------------------------------
% Created by : Celine Provins (05.2020)
%
% Load the GM regions in .mat format. If atlas exists only in .nii, save in .mat 
% separately each region of the atlas in the dimension dictated by cap_model.
%
% INPUT
% Savename is a string containing the prefix of the atlas for saving
%   purposes
% dirs is the struct extracted when loading the nifti file of the atlas using dir
% cap_model is a string containing the name of the model dictating the
%   dimension at which the masks extracted from the atlas are gonna be saved
%
% OUTPUT
% regions is a cell (1 x #regions). Each cell contains one region of the atlas
%   as a 1D vector.
%--------------------------------------------------------------------------
    maskPath = dirs.folder;
    if ~exist(fullfile(maskPath,strcat(savename,'_regions.mat')),'file')
        % Process the atlas and save the regions in .mat format
        disp(strcat('Processing atlas ',savename))
        
        %redimension atlas to match dimensions of cap_model
        cap_info = spm_vol(cap_model);
        dim = cap_info.dim;
        mask_name = dirs.name(1:end-4);
        mask = mapVolumeToVolume(strcat(mask_name,'.nii'),cap_model);
        %redimension atlas as a 1D vector
        mask = int16(mask(:));
        
        % create mask with only one region
        if strcmp('Yeo',savename)
            nbr_regions = 17; %18 to 34 are WM regions
        else
            nbr_regions = max(mask);
        end
        regions = cell(1,nbr_regions);
        for i = 1:nbr_regions
            idx_region = find(mask == i);
            %Each region is a 1D vector
            mask_region = zeros(1,dim(1)*dim(2)*dim(3));
            mask_region(idx_region)=1;
            regions{i}= mask_region;
        end
        save(fullfile(maskPath,...
                strcat(savename,'_regions.mat')), 'regions','-v7.3');
    else
        %Load the atlas
        region_struct = load(fullfile(maskPath,strcat(savename,'_regions.mat')));
        regions = region_struct.regions;
    end
end