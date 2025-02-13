function [active_PCC,idx_nonzero,idx_active] = active_PCC_unit(idx_frame,...
    isBC,prefix,suffix,defld_prefix,subject,BrainGraphPath,InpaintingPath,...
    subjectPath)
%--------------------------------------------------------------------------
% Modified by : Celine Provins (04.2020)
%
% Extract the 15% most active frames based on their activity in the PCC
%
% INPUT
% idx_frame : index of the most active frames per subject used only if isBC
%   is 0 because we retake for GSR the same frame indices that were
%   calculated for BC
% isBC : bool that precise if we work with BC data (1) or GSR (0)
% prefix : prefix of the file where the volumes are stored
% suffix : suffix of the file where the volumes are stored
% defld_prefix : prefix of the deformation field
% subject : name of the subject we are working with
% BrainGraphPath : path towards the folder containing the brain graphs
% InpaintingPath : path toward the folder containing the inpainted volumes
% subjectPath : path toward the folder containing the preprocessed volumes
%
% OUTPUT
% active_PCC : matrix containing the normalized active frames
% idx_nonzero : vector containing the linear coordinates of the non-zero
%   voxels in the normalized frames
% idx_active : vector containing the indices of the frames considered
%   active
%--------------------------------------------------------------------------
    %% Load variables
    indices_wb = load(fullfile(BrainGraphPath,subject,'ODF_Neigh_3_ODFPower_40',...
        'Spectrum_WB_improved','indices_wb.mat'));
    %non-zero indices in the inpainted volumes
    indices_wb = indices_wb.indices_wb;
    %Load the matrix containing the volumes (BC or GSR)
    subInpaintingPath = fullfile(InpaintingPath,subject);
    cd(subInpaintingPath);
    file = dir(strcat(prefix,'*',suffix));
    V = load(file.name);
    V = V.V;
    %matrix determining the dimension of the nifti to be saved
    fHeader = load(fullfile(InpaintingPath,subject,'fHeader.mat'));
    fHeader = fHeader.fHeader;

    %% Find the atlas normalized to the subject space
    segPath = fullfile(subjectPath,'struct','Segmented');
    cd(segPath);
    MyFolderInfo = dir('wNorm_AAL90_correctLR*');
    subject_atlas_path = MyFolderInfo.name;

    %% Detrend the volumes if it hasn't been done before
    if isBC
        disp('Let us detrend volume..')
        CUTNUMBER=10;
        SegmentLength = ceil(size(V,2) / CUTNUMBER);
        for iCut=1:CUTNUMBER
            if iCut~=CUTNUMBER
                Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
            else
                Segment = (iCut-1)*SegmentLength+1 : size(V,2);
            end
            V(:,Segment) = detrend(V(:,Segment));
            V(:,Segment) = detrend(V(:,Segment),'constant');
            fprintf('.');
        end
        V= zscore(V);
    end

    %% Reshape the brain graph to 4D volume
    epi = zeros([fHeader.dim,size(V,1)]);
    for i = 1:size(V,1)
        A=zeros(fHeader.dim);
        A(indices_wb) = V(i,:);
        epi(:,:,:,i) = A;
    end

    %% Extract the indices of the atlas corresponding to the PCC region
    mask = spm_read_vols(spm_vol(subject_atlas_path));
    idx_PCC_1 = find(mask == 35);
    idx_PCC_2 = find(mask == 36);
    idx_PCC = union(idx_PCC_1,idx_PCC_2);

    %% Extract maximum 15% active frame for PCC
    if isBC
        PCC_activity_mean = zeros(1,size(V,1));
        for i = 1:size(V,1)
            V_i = epi(:,:,:,i);
            PCC_activity_mean(i) = mean(V_i(idx_PCC));
        end
        perc_15_thres = prctile(PCC_activity_mean, 85);
        idx_active = find(PCC_activity_mean>=perc_15_thres);
    else
        %in case of GSR just take the same frame determined active in the BC case
        idx_active = idx_frame; 
    end

    %% Normalised the active frames to MNI space
    if exist(fullfile(BrainGraphPath,subject,'T1w','new_brainmask.nii'),'file')
        hdr_subject=cbiReadNiftiHeader(fullfile(BrainGraphPath,subject,'T1w','new_brainmask.nii'));
    else
        hdr_subject = cbiReadNiftiHeader(subject_atlas_path);
    end

    cd(segPath);
    defld = dir(defld_prefix);
    defld = fullfile(segPath,defld.name);
    [active_PCC,idx_nonzero] = NormMD(isBC,idx_active, V, defld,indices_wb,hdr_subject,fHeader,subInpaintingPath);
end