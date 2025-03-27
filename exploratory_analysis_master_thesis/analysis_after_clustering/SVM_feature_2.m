% Created by : Celine Provins (05.2020)
% Uses the files obtained and saved by SVM_feature_1.m
% The folder 'SVM-RFE-CBR-v1.3' containing the SVM-RFE code needs to be
% downloaded on https://ch.mathworks.com/matlabcentral/fileexchange/50701-feature-selection-with-svm-rfe
clear all; close all;
warning('off');

WholeBrain = 1; %work with WholeBrain CAPs (1) or PCC-seed CAPs (0)
numsel = 30; % 30 for resting state, 57 for working memory

%% Set up path
% celinePath = 'J:\Anjali_Diffusion_Pipeline\Celine';   
celinePath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
dataPath = fullfile(celinePath, 'data');
agccPath = fullfile(dataPath, 'RestingState');
ctrlPath = fullfile(dataPath, 'ControlsRS');
codePath = fullfile(celinePath,'analysis_after_clustering');
maskPathWM = fullfile(codePath,'Atlases','WMatlas');
maskPathGM = fullfile(codePath,'Atlases');
if WholeBrain
    capPathA = fullfile(agccPath, 'WholeBrainCAP');
    capPathC = fullfile(ctrlPath, 'WholeBrainCAP');
    capPathBoth = fullfile(dataPath,'WholeBrainCAP_RS');
    savePath = fullfile(codePath,'ResultsRFE_withACPC','WholeBrain');
    avgPath = fullfile(codePath,'ResultsClassification_withACPC','WholeBrain');
else
    capPathA = fullfile(agccPath, 'CAP');
    capPathC = fullfile(ctrlPath, 'CAP');
    capPathBoth = fullfile(dataPath,'CAP_RS_zscore');
    savePath = fullfile(codePath,'ResultsRFE_withACPC');
    avgPath = fullfile(codePath,'ResultsClassification_withACPC');
end
spmPath = fullfile(celinePath,'spm12');

addpath(genpath(spmPath));
addpath(genpath(codePath));

numsel = 30; %nbr of frame selected per subject

%% Load variables
% model that dictates the dimension of the files to be saved
cap_model = fullfile(capPathBoth,'GSR_CAP_1_mean.nii');
cap_info = spm_vol(cap_model);
dim = cap_info.dim;

accuracy_struct = load(fullfile(savePath, 'accuracy_opt_default.mat'));
accuracySVM_WM = accuracy_struct.mean_acc_WM;
accuracySVM_GM = accuracy_struct.mean_acc_GM;

%Pairs to which CAP each frame contributes
frame_check = load(fullfile(capPathBoth,'frame_check.mat'));
frame = frame_check.frame_GSR;

%Active frames from AgCC
activeA = load(fullfile(capPathA,'active.mat'));
activeA_truncated=activeA.active_GSR_truncated;
idxA = activeA.idx;

%Active frames from controls
activeC = load(fullfile(capPathC,'active.mat'));
activeC_truncated=activeC.active_GSR_truncated;
idxC = activeC.idx;
idx = union(idxA, idxC);

%define constants
NC = size(activeC_truncated,1); %nbr of frames from controls
NA = size(activeA_truncated,1); %nbr of frames from AgCC
nframes = NA + NC;

%Build the vector that indicates to which subject belong the frame i
numsel_vec = build_numsel_vec(WholeBrain,activeA,activeC,NC,nframes,numsel);

%Assemble the active frames from controls and AgCC
active = zeros(nframes,dim(1)*dim(2)*dim(3));
for i = 1:NC
    active(i,idxC) = activeC_truncated(i,:);
end
for i = 1:NA
    active(i+NC,idxA) = activeA_truncated(i,:);
end
active = active(:,idx);
active = zscore(active,0,2);

%Find nbr of clusters
delimiter = find(frame==0); 
%the indices corresponding to a cluster are separated by zeros in the vector
k = length(delimiter)+1; %nbr of clusters+1(all)  

%% Load the average and masked frames
if WholeBrain
    avg_struct = load(fullfile(avgPath,strcat('average_within_mask_WB.mat')));
    avg_GM = avg_struct.avg_YEO;
else
    avg_struct = load(fullfile(avgPath,strcat('average_within_mask.mat')));
    avg_GM = avg_struct.avg_AH;
end
avg_WM = avg_struct.avg_WM;
nbr_regionGM = size(avg_GM,2);
nbr_regionWM = size(avg_WM,2);

%% Extracting the frames indices contributing to each CAP
list_idx_cap_nonordered = {1,k}; %necessary because CAP_RS and CAP_RS_zscore not same #CAP
list_idx_cap = {1,k};
for l = 1:k
    %Identify the index of frames contributing to CAP-l
    if l==k
        list_idx_cap_nonordered{l} = 1:nframes; %take all frames  
    elseif l==(k-1)
        list_idx_cap_nonordered{l} = frame(delimiter(l)+1:end); 
    else
        list_idx_cap_nonordered{l} = frame(delimiter(l)+1:delimiter(l+1)-1);
    end
end
%correspondance between #CAP of zscored to non-zscored (#CAP on visualization image)
list_idx_cap{1} = list_idx_cap_nonordered{1};
list_idx_cap{2} = list_idx_cap_nonordered{6};
list_idx_cap{3} = list_idx_cap_nonordered{7};
list_idx_cap{4} = list_idx_cap_nonordered{3};
list_idx_cap{5} = list_idx_cap_nonordered{5};
list_idx_cap{6} = list_idx_cap_nonordered{8};
list_idx_cap{7} = list_idx_cap_nonordered{2};
list_idx_cap{8} = list_idx_cap_nonordered{4};
list_idx_cap{9} = list_idx_cap_nonordered{9};
    
%% Feature extraction
ftRankWM = zeros(k,nbr_regionWM);
ftRankGM = zeros(k,nbr_regionGM);
mean_acc_WM = zeros(1,k);
mean_acc_GM = zeros(1,k);

threshold = [1,0.9,0.8,0.7,0.6,0.5];
for l = 1:k 
    disp(strcat('Analyis for CAP',int2str(l)))

    %Identify the index of frames contributing to CAP-l
    idx_cap = list_idx_cap{l};

    %Restrict the matrices to the frames contributing to cap l
    avg_WM_truncated = avg_WM(idx_cap,:);
    avg_GM_truncated = avg_GM(idx_cap,:);

    %Find which subjects contribute to that cap
    subjects_cap = [];
    for i = 1:length(idx_cap)
        m = numsel_vec(idx_cap(i));
        subjects_cap = [subjects_cap m];
    end  
    %so that each subject appear only once in the list
    subjects_cap = unique(subjects_cap);
    labels = label(idx_cap,NC);
%     save(fullfile(savePath,'labels.mat'),'labels');
%     save(fullfile(savePath,'avg_WM.mat'),'avg_WM_truncated');

    %select threshold with best accuracy
    [~,tWM] = max(accuracySVM_WM(:,l));
    [~,tGM] = max(accuracySVM_GM(:,l));
    t=0;
    %Call SVM-RFE
    [ftRankWM(l,:),~] = ftSel_SVMRFECBR_ori(zscore(avg_WM_truncated),labels',l,t,threshold(tWM),savePath,1,WholeBrain);
    [ftRankGM(l,:),~] = ftSel_SVMRFECBR_ori(zscore(avg_GM_truncated),labels',l,t,threshold(tGM),savePath,0,WholeBrain);
end

save(fullfile(savePath,'RFE_results_default.mat'),'ftRankGM','ftRankWM');