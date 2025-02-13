function idx_active = select_active_frames(Vnorm,subjectPath,nii_model)
%--------------------------------------------------------------------------
% Created by : Celine Provins (06.2020)
%
% For each of the region of the YEO atlas, select the 15% most active
% frames based on their activity in this particular region. The final
% selection of active frames corresponds to the union of the most active
% frames of each region.
%
% INPUT
% Vnorm :
% subjectPath : path towards the preprocessed data of a particular subject
% nii_model : foldername + filename of the model nifti that determines the 
%   dimensions at which the atlas is resampled
%
% OUTPUT
% idx_active : vector containing the indices of the frames we consider active 
%--------------------------------------------------------------------------

    %% Load atlas
    segPath = fullfile(subjectPath,'struct','Segmented');
    cd(segPath);
    atlas_struct = dir('wNorm_WM_atlas_Yeo_Combined17*.nii');
    subject_atlas_path = atlas_struct.name;
    %Resample the atlas so that its dimensions match with the dimensions of
    %the model
    atlas = mapVolumeToVolume(subject_atlas_path,nii_model);
    
    nbr_region = 17;

    %% Extract 15% most active frames based on activity of each region
    idx_active = [];
    for r = 1:nbr_region
        %Find indices corresponding to one region
        idx_region = find(atlas == r);

        % extract maximum 15% active frame for the region
        nbr_frames = size(Vnorm,1);
        region_activity_mean = zeros(1,nbr_frames);
        for i = 1:nbr_frames
            region_activity_mean(i) = mean(Vnorm(i,idx_region));
        end
        perc_15_thres = prctile(region_activity_mean, 85);
        idx_active = [idx_active find(region_activity_mean>=perc_15_thres)];   	
    end
    %keep only once all the indices
    idx_active = unique(idx_active);
end