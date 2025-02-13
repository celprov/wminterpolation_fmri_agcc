function whole_brain_cap_one(dataset)
%--------------------------------------------------------------------------
% Created by : Celine Provins (06.2020)
%
% Extract the whole-brain CAPs for the control or the AgCC condition separated
%
% dataset : string that precise with which data you want to work
%--------------------------------------------------------------------------

    %% Set up the paths
%     dataset = 'RestingState';
    warning('off');
    % celinePath = 'J:\Anjali_Diffusion_Pipeline\Celine';   
    celinePath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
    BrainGraphPath = '/media/miplab-nas2/Data/AgCC/Anjali/BrainGraph_results';
    spmPath = fullfile(celinePath,'spm12');
    codeBasePath = fullfile(celinePath,'WholeBrainCaps');
    codeBasePath2 = fullfile(celinePath,'kmeans');
    cbiPath = fullfile(celinePath,'mrTools');
    dataPath = fullfile(celinePath,'data',dataset);
    InpaintingPath = fullfile(dataPath,'Inpainting_results');
    capPath = fullfile(dataPath,'WholeBrainCAP');
    if ~exist(capPath,'dir')
        mkdir(dataPath,'WholeBrainCAP');
    end

    addpath(genpath(spmPath));
    addpath(genpath(codeBasePath));
    addpath(genpath(codeBasePath2));
    addpath(genpath(cbiPath));
    
    %% Selecting patient folders
    dirs = dir(fullfile(dataPath, 's*'));
    subjects = cell(size(dirs));
    for i = 1:size(dirs,1)
        subjects{i} = dirs(i).name;
    end
    
    %% Defining parameters
    deformationPrefix = 'y_*';
    idx_active = cell(1,length(subjects));
    nii_model = fullfile(InpaintingPath,subjects{1},'w3WholeBrainBC001.nii');
    
    %% BC Volume
    isBC = 1;
    prefix = 'BCVolumes_Cov_';
    active_BC = [];
    active_BC_lr = []; %low_resolution volumes to reduce computation time of kmean

    for i = 1:length(subjects)
        fprintf('%s%d%s%d\n','BC : Processing subject number ', i,...
            '. Total number of subjects is ', size(subjects,1));
        subject = subjects{i};
        subjectPath = fullfile(dataPath,subject); 
        segPath = fullfile(subjectPath,'struct','Segmented');
        
        %Normalize to MNI and resize to lower resolution ALL the volumes
        [Vnorm,VnormLR,idx_nonzeroBC,idx_nonzeroBC_lr] = Norm_InpaintVol(isBC, prefix,...
            deformationPrefix,BrainGraphPath,InpaintingPath,subject,segPath);
        
        idx_active_sub = select_active_frames(VnormLR,subjectPath,nii_model);
        idx_active{i} = idx_active_sub;
        
        active_BC = vertcat(active_BC,Vnorm(idx_active_sub,:));
        active_BC_lr = vertcat(active_BC_lr,VnormLR(idx_active_sub,:));
    end
%     threshold = 0.7;
%     idx_BC = select_indices(active_BC,threshold);
    
    clearvars Vnorm
    clearvars VnormLR
    
    %% GSR volume
    isBC = 0;
    prefix = 'GSR_Lambda';
    active_GSR = [];
    active_GSR_lr = []; %low_resolution volumes to reduce computation time of kmean
    
    for i = 1:length(subjects)
        fprintf('%s%d%s%d\n','GSR : Processing subject number ', i,...
            '. Total number of subjects is ', size(subjects,1));
        subject = subjects{i};
        subjectPath = fullfile(dataPath,subject); 
        segPath = fullfile(subjectPath,'struct','Segmented');
        
        %Normalize to MNI and resize to lower resolution ALL the inpainted volumes
        [Vnorm,VnormLR,idx_nonzeroGSR,idx_nonzeroGSR_lr] = Norm_InpaintVol(isBC, prefix,...
            deformationPrefix,BrainGraphPath,InpaintingPath,subject,segPath);
        
        active_GSR = vertcat(active_GSR, Vnorm(idx_active{i},:));
        active_GSR_lr = vertcat(active_GSR_lr, VnormLR(idx_active{i},:));
    end
%     threshold = 0.7;
%     idx_GSR = select_indices(active_GSR,threshold);
    
    clearvars Vnorm
    clearvars VnormLR
    
    %% k_means 
    %Truncate matrices + cluster in low resolution to improve speed
    idx_lr = union(idx_nonzeroGSR_lr,idx_nonzeroBC_lr);
    active_GSR_truncated_lr = active_GSR_lr(:,idx_lr);
    active_BC_truncated_lr = active_BC_lr(:,idx_lr);
    
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
    cluster_idx_BC = kmeans(active_BC_truncated_lr,k,'replicate',10,'Distance','cosine');
    cluster_idx_GSR = kmeans(active_GSR_truncated_lr,k,'replicate',10,'Distance','cosine');
    
    save(fullfile(capPath,'active_lowResolution.mat'),'idx_lr',...
        'active_BC_truncated_lr','active_GSR_truncated_lr','-v7.3');
    
    clearvars active_BC_truncated_lr
    clearvars active_GSR_truncated_lr
    clearvars idx_lr

    %% align and save the CAPs
    %Truncate matrices + obtain CAPs in high resolution
    idx = union(idx_nonzeroGSR,idx_nonzeroBC);
    active_GSR_truncated = active_GSR(:,idx);
    active_BC_truncated = active_BC(:,idx);
    
    disp(strcat('Getting the CAPs for ', dataset));
    hdr = cbiReadNiftiHeader(fullfile(InpaintingPath,subject,'w2WholeBrainIN001.nii'));
    CAP_GSR = getCAPs(active_GSR_truncated,k,cluster_idx_GSR);
    CAP_BC = getCAPs(active_BC_truncated,k,cluster_idx_BC);
    [CAP_GSR,idx_col] = arange_CAPs(CAP_BC,CAP_GSR);

    disp(strcat('Saving the CAPs for ',dataset));
    saveCAPs('BC_',CAP_BC,hdr,k,idx,1,idx_col,capPath);
    saveCAPs('GSR_',CAP_GSR,hdr,k,idx,0,idx_col,capPath);

    %% Save frames in group wrt to their CAP membership
    frame_BC = [];
    frame_GSR = [];
    for i = 1:k
        frame_number_BC = find(cluster_idx_BC==i);
        frame_number_GSR = find(cluster_idx_GSR==idx_col(i));
        frame_BC = vertcat(frame_BC,0,frame_number_BC);
        frame_GSR = vertcat(frame_GSR,0,frame_number_GSR);
    end
    save(fullfile(capPath,'frame_check.mat'),'frame_BC','frame_GSR','idx_active');

    %% save variables
    save(fullfile(capPath,'active.mat'),'idx','idx_active','active_BC_truncated','active_GSR_truncated','-v7.3')
    save(fullfile(capPath,'CAP.mat'),'idx','CAP_BC','CAP_GSR','-v7.3')
end

