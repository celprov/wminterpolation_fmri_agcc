% function extractCAP_both(RestingState_mode)
%--------------------------------------------------------------------------
% Modified by : Celine Provins (04.2020)
%
% Extract the PCC-seed CAPs using both control and AgCC subjects combined
%
% RestingState_mode : bool precising whether you want to work with resting
%   state data (1) or working memory data (0)
%--------------------------------------------------------------------------
    RestingState_mode = 1;
    
    %% Defining the paths
    %celinePath = 'J:\Anjali_Diffusion_Pipeline\Celine';
    celinePath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
    spmPath = fullfile(celinePath,'spm12');
    codeBasePath = fullfile(celinePath, 'kmeans');
    cbiPath = fullfile(celinePath,'mrTools');
    dataPath = fullfile(celinePath, 'data');
    if RestingState_mode
        InpaintingPath = fullfile(dataPath,'RestingState','Inpainting_results');
        capPath_AgCC = fullfile(dataPath,'RestingState','CAP_test_RemoveMean');
        capPath_Ctrl = fullfile(dataPath,'ControlsRS','CAP_test_RemoveMean');
        BothCapPath = fullfile(dataPath,'CAP_RS_test_zscore');
        if ~exist(BothCapPath,'dir')
        mkdir(dataPath,'CAP_RS_test_zscore');
        end
    else
        InpaintingPath = fullfile(dataPath,'WorkingMemory','Inpainting_results');
        capPath_AgCC = fullfile(dataPath,'WorkingMemory','CAP');
        capPath_Ctrl = fullfile(dataPath,'ControlsWM','CAP');
        BothCapPath = fullfile(dataPath,'CAP_WM');
        if ~exist(BothCapPath,'dir')
        mkdir(dataPath,'CAP_WM');
        end
    end
    
    addpath(genpath(spmPath));
    addpath(genpath(codeBasePath));
    addpath(genpath(cbiPath));
    %% load the data generated with AgCC and Ctrl separated

    active_AgCC = load(fullfile(capPath_AgCC,'active.mat'));
    active_AgCC_BC_truncated = active_AgCC.active_BC_truncated;
    active_AgCC_GSR_truncated = active_AgCC.active_GSR_truncated;
    idx_AgCC = active_AgCC.idx;

    active_Ctrl = load(fullfile(capPath_Ctrl,'active.mat'));
    active_Ctrl_BC_truncated = active_Ctrl.active_BC_truncated;
    active_Ctrl_GSR_truncated = active_Ctrl.active_GSR_truncated;
    idx_Ctrl = active_Ctrl.idx;

    idx = union(idx_Ctrl,idx_AgCC);
    NC = size(active_Ctrl_BC_truncated,1); %nbr of control frames
    NA = size(active_AgCC_BC_truncated,1); %nbr of AgCC frames
    N = NA+NC;
    % Nifti model that determines dimension of saved CAPs
    if RestingState_mode
        hdr=cbiReadNiftiHeader(fullfile(InpaintingPath,'s012','wINsel001.nii'));
    else
        hdr=cbiReadNiftiHeader(fullfile(InpaintingPath,'s102','wINsel001.nii'));
    end
    %% BC active frames
    V_BC = zeros(N,hdr.dim(2)*hdr.dim(3)*hdr.dim(4));
    disp('Extracting the BC data');
    for i = 1:NC
        V_BC(i,idx_Ctrl) = active_Ctrl_BC_truncated(i,:);
    end
    for i = 1:NA
        V_BC(i+NC,idx_AgCC) = active_AgCC_BC_truncated(i,:);
    end
    V_BC = V_BC(:,idx);

    %% GSR active frames
    V_GSR = zeros(N,hdr.dim(2)*hdr.dim(3)*hdr.dim(4));
    disp('Extracting the GSR data');
    for i = 1:NC
        V_GSR(i,idx_Ctrl) = active_Ctrl_GSR_truncated(i,:);
    end
    for i = 1:NA
        V_GSR(i+NC,idx_AgCC) = active_AgCC_GSR_truncated(i,:);
    end
    V_GSR = V_GSR(:,idx);

    %% K-MEANS
    optimize_k = 0; %% decide whether to select an optimal k
    k_range = [2,3,4,5,6,7,8];
    if optimize_k
        disp('Finding optimal nbr of clusters');
        [k_opt] = find_optimal_k(active_Ctrl_BC_truncated,k_range,BothCapPath);
        k = k_opt;
    else
        k = 8;
    end
    disp('Performing k-means');
    cluster_idx_BC = kmeans(V_BC,k,'replicate',10,'Distance','cosine');
    cluster_idx_GSR = kmeans(V_GSR,k,'replicate',10,'Distance','cosine');

    %% align and save the CAPs

    disp('Getting the CAPs for both AgCC and Controls combined');
    CAP_both_BC = getCAPs(V_BC,k,cluster_idx_BC);
    CAP_both_GSR = getCAPs(V_GSR,k,cluster_idx_GSR);
    %Match the number of CAP_GSR with CAP_BC
    [CAP_both_GSR,idx_GSR] = arange_CAPs(CAP_both_BC,CAP_both_GSR);

    disp('Saving the CAPs for both AgCC and Controls combined');
    saveCAPs('BC_',CAP_both_BC,hdr,k,idx,1,idx_GSR,BothCapPath);
    saveCAPs('GSR_',CAP_both_GSR,hdr,k,idx,0,idx_GSR,BothCapPath);

    %% Save frames in group wrt to their CAP membership
    frame_BC = [];
    frame_GSR = [];
    for i = 1:k
        frame_BC = vertcat(frame_BC,0,find(cluster_idx_BC==i));
        frame_GSR = vertcat(frame_GSR,0,find(cluster_idx_GSR==idx_GSR(i)));
    end
    save(fullfile(BothCapPath,'frame_check.mat'),'frame_BC','frame_GSR');
    save(fullfile(BothCapPath,'CAP.mat'),'idx','CAP_both_BC','CAP_both_GSR','-v7.3');%,'WM_MASK','mask_GM','-v7.3')

    %% find frames distribution in each CAP
    V = zeros(hdr.dim(2),hdr.dim(3),hdr.dim(4));
    BC_dist = zeros(1,k);
    GSR_dist = zeros(1,k);
    for i = 1:k
        %find all indices that belongs to cluster i
        idx_BC_cluster_i = find(cluster_idx_BC==i);
        idx_GSR_cluster_i = find(cluster_idx_GSR==i);
        %find all indices from control that belongs to cluster i
        BC_dist(i) = size(find(idx_BC_cluster_i<NC))/size(idx_BC_cluster_i);
        GSR_dist(i) = size(find(idx_GSR_cluster_i<NC))/size(idx_GSR_cluster_i);
    end
    %find all indices from AgCC that belongs to cluster i
    BC_dist_AgCC = ones(1,k) - BC_dist;
    GSR_dist_AgCC = ones(1,k) - GSR_dist;
    save(fullfile(BothCapPath,'distribution_data.mat'),'BC_dist_AgCC',...
        'GSR_dist_AgCC','BC_dist','GSR_dist');
% end