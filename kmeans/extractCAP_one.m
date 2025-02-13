function extractCAP_one(dataset, num_frame_sel)
%--------------------------------------------------------------------------
% Modified by : Celine Provins (04.2020)
%
% Extract the PCC-seed CAPs for the control or the AgCC condition separated
%
% dataset : string that precise with which data you want to work
% num_frame_sel : integer that precise how many active frames per subject
%   where selected
%--------------------------------------------------------------------------
%     dataset = 'RestingState';
%     num_frame_sel =30;
    warning('off');
    
    %% Defining the paths
    % celinePath = 'J:\Anjali_Diffusion_Pipeline\Celine';   
    celinePath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
    spmPath = fullfile(celinePath,'spm12');
    codeBasePath = fullfile(celinePath,'kmeans');
    cbiPath = fullfile(celinePath,'mrTools');
    BrainGraphPath = '/media/miplab-nas2/Data/AgCC/Anjali/BrainGraph_results';
    dataPath = fullfile(celinePath,'data',dataset);
    InpaintingPath = fullfile(dataPath,'Inpainting_results');
    capPath = fullfile(dataPath,'CAP');
    if ~exist(capPath,'dir')
        mkdir(dataPath,'CAP');
    end

    addpath(genpath(spmPath));
    addpath(genpath(codeBasePath));
    addpath(genpath(cbiPath));
    
    %% Selecting patient folders
    dirs = dir(fullfile(dataPath, 's*'));
    subjects = cell(size(dirs));
    for i = 1:size(dirs,1)
        subjects{i} = dirs(i).name;
    end
    
    %% Defining parameters
    deformationPrefix = 'y_*'; %prefix of the deformation file
    suffix = '.mat';
    idx_active= zeros(length(subjects) * num_frame_sel,1);  % index of the most active frames per subject: column - subjects
    
    %% BC Volume : CAP using grey matter only 
    isBC = 1;
    prefix = 'BCVolumes_Cov_';
    active_BC=[];
    for i = 1:length(subjects)
        fprintf('%s%d%s%d\n','BC : Processing subject number ', i,...
            '. Total number of subjects is ', size(subjects,1));
        subject = subjects{i};
        [active,idx_nonzero,idx_act] = active_PCC_unit(...
            idx_active(1+num_frame_sel*(i-1):num_frame_sel*i,:),...
            isBC,prefix,suffix,deformationPrefix,subject,BrainGraphPath,...
            InpaintingPath,fullfile(dataPath,subject));
        idx_active(1+num_frame_sel*(i-1):num_frame_sel*i,:) = idx_act;
        active_BC = vertcat(active_BC,active);
    end
    %Since volumes are not perfecty aligned, keep voxels that appear
    %non-zero in at least 70% of the frames
    threshold = 0.7;
    idx_BC = select_indices(active_BC,threshold);
    
    clearvars active
    clearvars idx_nonzero
    clearvars idx_act  
    save(fullfile(capPath,'active_BC.mat'),'active_BC','idx_BC','-v7.3');
    %% GSR volume : CAP using inpainted volume
    isBC = 0;
    prefix = 'GSR_Lambda';
    active_GSR=[];
    for i = 1:length(subjects)
        fprintf('%s%d%s%d\n','GSR : Processing subject number ', i,...
            '. Total number of subjects is ', size(subjects,1));
        subject = subjects{i};
        [active,~] = SCA_active_PCC_unit(...
            idx_active(1+num_frame_sel*(i-1):num_frame_sel*i,:),...
            isSCA,isBC,prefix,suffix,deformationPrefix,subject,BrainGraphPath,...
            InpaintingPath,fullfile(dataPath,subject));
        active_GSR = vertcat(active_GSR,active);
    end
    %Since volumes are not perfecty aligned, keep voxels that appear
    %non-zero in at least 70% of the frames
    threshold = 0.7;
    idx_GSR = select_indices(active_GSR,threshold);
    clearvars active;
    save(fullfile(capPath,'active_GSR.mat'),'active_GSR','idx_GSR','-v7.3');
    
    %% Truncate data
    idx = union(idx_GSR,idx_BC);
    active_GSR_truncated = active_GSR(:,idx);
    active_BC_truncated = active_BC(:,idx);

    %% k_means 
    optimize_k = 0; %% decide whether to select an optimal k
    k_range = [2,3,4,5,6,7,8];
    if optimize_k
        disp('Find optimal nbr of clusters');
        k_opt = find_optimal_k(active_BC_truncated,k_range,capPath);
        k = k_opt;
    else
        k = 8;
    end
    disp('Executing K-means');
    cluster_idx_BC = kmeans(active_BC_truncated,k,'replicate',10,'Distance','cosine');
    cluster_idx_GSR = kmeans(active_GSR_truncated,k,'replicate',10,'Distance','cosine');

    %% align and save the CAPs
    disp(strcat('Getting the CAPs for ', dataset));
    %Nifti model that determines the dimension of the CAPs
    hdr=cbiReadNiftiHeader(fullfile(InpaintingPath,subject,'wINsel001.nii'));
    CAP_GSR = getCAPs(active_GSR_truncated,k,cluster_idx_GSR);
    CAP_BC = getCAPs(active_BC_truncated,k,cluster_idx_BC);
    %Match the number of CAP_GSR with CAP_BC
    [CAP_GSR,idx_GSR] = arange_CAPs(CAP_BC,CAP_GSR);

    disp(strcat('Saving the CAPs for ',dataset));
    saveCAPs('BC_',CAP_BC,hdr,k,idx,1,idx_GSR,capPath);
    saveCAPs('GSR_',CAP_GSR,hdr,k,idx,0,idx_GSR,capPath);

    %% Save frames in group wrt to their CAP membership
    frame_BC = [];
    frame_GSR = [];
    for i = 1:k
        frame_number_BC = find(cluster_idx_BC==i);
        frame_number_GSR = find(cluster_idx_GSR==idx_GSR(i));
        frame_BC = vertcat(frame_BC,0,frame_number_BC);
        frame_GSR = vertcat(frame_GSR,0,frame_number_GSR);
    end
    save(fullfile(capPath,'frame_check.mat'),'frame_BC','frame_GSR','idx_active');

    %% save variables
    save(fullfile(capPath,'active.mat'),'idx','active_BC_truncated','active_GSR_truncated','-v7.3')
    save(fullfile(capPath,'CAP.mat'),'idx','CAP_BC','CAP_GSR','-v7.3')
end

