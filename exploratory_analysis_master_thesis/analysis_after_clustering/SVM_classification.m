% Created by : Celine Provins (05.2020)
clear all; close all;
warning('off');

WholeBrain = 1; %work with WholeBrain CAPs (1) or PCC-seed CAPs (0)
numsel = 30; % 30 for resting state, 57 for working memory

%% Set up paths
%celinePath = 'J:\Anjali_Diffusion_Pipeline\Celine';   
celinePath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
dataPath = fullfile(celinePath, 'data');
agccPath = fullfile(dataPath, 'RestingState');
ctrlPath = fullfile(dataPath, 'ControlsRS');
codePath = fullfile(celinePath,'analysis_after_clustering');
codePath2 = fullfile(celinePath,'kmeans');
maskPathWM = fullfile(codePath,'Atlases','WMatlas');
maskPathGM = fullfile(codePath,'Atlases');
cbiPath = fullfile(celinePath,'mrTools');
spmPath = fullfile(celinePath,'spm12');
if WholeBrain
    capPathA = fullfile(agccPath, 'WholeBrainCAP');
    capPathC = fullfile(ctrlPath, 'WholeBrainCAP');
    capPathBoth = fullfile(dataPath,'WholeBrainCAP_RS');
    savePath = fullfile(codePath, 'ResultsClassification_withACPC','WholeBrain');
    RFEPath = fullfile(codePath, 'ResultsRFE_withACPC','WholeBrain');
else
    capPathA = fullfile(agccPath, 'CAP');
    capPathC = fullfile(ctrlPath, 'CAP');
    capPathBoth = fullfile(dataPath,'CAP_RS_zscore');
    savePath = fullfile(codePath, 'ResultsClassification_withACPC');
    RFEPath = fullfile(codePath, 'ResultsRFE_withACPC');
end
addpath(genpath(spmPath));
addpath(genpath(cbiPath));
addpath(genpath(codePath));
addpath(genpath(codePath2));

%% Load variables
%Model dictating the dimension of the files to be saved
cap_model = fullfile(capPathBoth,'GSR_CAP_1_mean.nii');
cap_info = spm_vol(cap_model);
dim = cap_info.dim;

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

% cap_struct = load(fullfile(capPathBoth,'CAP.mat'));
% cap = cap_struct.CAP_both_GSR;
% indices_cap = cap_struct.idx;

%define constants
NC = size(activeC_truncated,1); %nbr of frames from controls
NA = size(activeA_truncated,1); %nbr of frames from AgCC
nframes = NA + NC;

%Build the vector that indicates to which subject belong the frame i
numsel_vec = build_numsel_vec(WholeBrain,activeA,activeC,NC,nframes,numsel);

%Assemble the active frames from control and AgCC
active = zeros(nframes,dim(1)*dim(2)*dim(3));
for i = 1:NC
    active(i,idxC) = activeC_truncated(i,:);
end
for i = 1:NA
    active(i+NC,idxA) = activeA_truncated(i,:);
end
active = active(:,idx);
active = zscore(active,0,2);

%% Finding each patient folder
dirsA = dir(fullfile(agccPath, 's*'));
dirsC = dir(fullfile(ctrlPath, 's*'));
na = size(dirsA,1); %number of AgCC subjects
nc = size(dirsC,1); %number of Control subjects
n= na + nc;
subjects = cell(1,n);
for i = 1:nc
    subjects{i} = fullfile(dirsC(i).folder,dirsC(i).name);
end
for i = 1:na
    subjects{i+nc} = fullfile(dirsA(i).folder,dirsA(i).name);
end

%% Plot the mask for nice visualization of CAPs
% InpaintingPath = fullfile(dataPath,'RestingState','Inpainting_results');
% subject = subjects{1};
% hdr=cbiReadNiftiHeader(fullfile(InpaintingPath,subject,'wINsel001.nii'));
% Vmask = ones(hdr.dim(2),hdr.dim(3),hdr.dim(4));
% Vmask(idx) = 0;
% cbiWriteNifti(fullfile(savePath,strcat('Vmask.nii')),Vmask,hdr,'float32');

%% Processing the atlases
cd(maskPathGM)
AH_dirs = dir('LH-Andrews-Hannah-atlas*.nii');
Yeo_dirs = dir('*Yeo*.nii');
AHregions = load_GM_atlas('AH', AH_dirs, cap_model);
YEOregions = load_GM_atlas('Yeo', Yeo_dirs, cap_model);
WMregions = load_WM_atlas(cap_model,maskPathWM,RFEPath);

%% Averaging the signal in the WM and GM mask

if ~exist(fullfile(savePath,strcat('average_within_maskWB.mat')),'file')
    %% Processing the WM and GM masks
    c2_mask = cell(1,length(subjects));
    c1_mask = cell(1,length(subjects));
    for i = 1:length(subjects)
        segPath = fullfile(subjects{i},'struct','Segmented');
        if ~exist(fullfile(segPath,'GM_mask.mat'),'file')
            %Processing the masks and save them in .mat format
            disp('Processing the masks')
            cd(segPath)
            c2str=dir('MNI_c2*');
            c1str=dir('MNI_c1*');
            %redimension the mask to match dimensions of CAPs
            c2 = mapVolumeToVolume(fullfile(c2str.folder,c2str.name),cap_model);
            c1 = mapVolumeToVolume(fullfile(c1str.folder,c1str.name),cap_model);
            %reshape to 1D
            c2_mask{i} = c2(:);
            c1_mask{i} = c1(:);

            save(fullfile(segPath,'WM_mask.mat'),'c2_mask');
            save(fullfile(segPath,'GM_mask.mat'),'c1_mask');
        else
            %Load the masks
            c2_struct = load(fullfile(segPath,'WM_mask.mat'));
            c2_mask = c2_struct.c2_mask;
            c1_struct = load(fullfile(segPath,'GM_mask.mat'));
            c1_mask = c1_struct.c1_mask;     
        end
    end  
    
    %% Computing average of each frame within the atlas regions
    %% Apply the WM and GM masks to the frames
    disp('Computing average within mask')
    avg_WM = zeros(nframes,length(WMregions));
    avg_GM = zeros(nframes,length(AHregions));
    avg_GM = zeros(nframes,length(YEOregions));
    active_maskedWM = zeros(size(active));
    active_maskedGM = zeros(size(active));

    for i = 1:nframes
        avg_WM(i,:) = average_within_mask(WMregions,idx,active(i,:));
        avg_AH(i,:) = average_within_mask(AHregions,idx,active(i,:));
        avg_YEO(i,:) = average_within_mask(YEOregions,idx,active(i,:));
        
        m = numsel_vec(i);
        WM_mask = c2_mask{m};
        GM_mask = c1_mask{m};
        active_maskedWM(i,:) = apply_mask(WM_mask,idx,active(i,:));
        active_maskedGM(i,:) = apply_mask(GM_mask,idx,active(i,:));
    end

    if WholeBrain
        save(fullfile(savePath,strcat('average_within_mask_WB.mat')),...
            'avg_WM', 'avg_AH','avg_YEO');
        save(fullfile(savePath,strcat('active_masked_WB.mat')),...
            'active_maskedWM','active_maskedGM','-v7.3');
    else
        save(fullfile(savePath,strcat('average_within_mask.mat')),...
            'avg_WM', 'avg_AH','avg_YEO');
        save(fullfile(savePath,strcat('active_masked.mat')),...
            'active_maskedWM','active_maskedGM','-v7.3');
    end
end

%% Load the average and the masked frames
if WholeBrain
    avg_struct = load(fullfile(savePath,strcat('average_within_mask_WB.mat')));
    avg_GM = avg_struct.avg_YEO;
    mask_struct = load(fullfile(savePath,strcat('active_masked_WB.mat')));

else
    avg_struct = load(fullfile(savePath,strcat('average_within_mask.mat')));
    avg_GM = avg_struct.avg_AH;
    mask_struct = load(fullfile(savePath,strcat('active_masked.mat')));
end
avg_WM = avg_struct.avg_WM;
active_maskedWM = mask_struct.active_maskedWM;
active_maskedGM = mask_struct.active_maskedGM;


%% Averaging the signal within the CAPs directly
% delimiter = find(frame==0);
% k = length(delimiter)+1; %nbr of clusters+1(all frames)
% %the list of indices specific to each cap are separated by zeroes
% avg_capWM = zeros(k-1,length(WMregions));
% avg_capYEO = zeros(k-1,length(YEOregions));
% for l = 1:k-1
%     avg_capWM(l,:) = average_within_mask(WMregions,indices_cap,cap(l,:));
%     avg_capYEO(l,:) = average_within_mask(YEOregions,indices_cap,cap(l,:));
% end

%% Extracting the frames indices contributing to each CAP
delimiter = find(frame==0);
%the list of indices specific to each cap are separated by zeroes
k = length(delimiter)+1; %nbr of clusters+1(all frames)

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

%% Find list of frames contributing to each cap
mean_acc = zeros(6,k);
avg_WM_both = zeros(k-1,size(avg_WM,2));
avg_WM_ctrl = zeros(k-1,size(avg_WM,2));
avg_WM_agcc = zeros(k-1,size(avg_WM,2));
avg_GM_both = zeros(k-1,size(avg_GM,2));
avg_GM_ctrl = zeros(k-1,size(avg_GM,2));
avg_GM_agcc = zeros(k-1,size(avg_GM,2));

confusion_mat_IN_sum = zeros(k,2,2);
confusion_mat_WM_sum = zeros(k,2,2);
confusion_mat_GM_sum = zeros(k,2,2);
confusion_mat_WMpca_sum = zeros(k,2,2);
confusion_mat_GMpca_sum = zeros(k,2,2);

WM_ctrl_list = cell(1,k-1);
WM_agcc_list = cell(1,k-1);
GM_ctrl_list = cell(1,k-1);
GM_agcc_list = cell(1,k-1);
        
for l= 1:k
    disp(strcat('Analyis for CAP',int2str(l)))
    
    %Identify the index of frames contributing to CAP-l
    idx_cap = list_idx_cap{l};
    nframes_cap = length(idx_cap);
    
    %Restrict the matrices to the frames contributing to cap l
    avg_WM_truncated = avg_WM(idx_cap,:);
    avg_GM_truncated = avg_GM(idx_cap,:);

    %% Find subjects contributing to each CAP
    subjects_cap = [];
    for i = 1:nframes_cap
            %Find which subjects contribute to that cap
            m = numsel_vec(idx_cap(i));
            subjects_cap = [subjects_cap m];
    end
    %so that each subject appear only once in the list
    subjects_cap = unique(subjects_cap);
    
    %% Extract average within atlas to plot heatmap
    if l~=k
        WM_ctrl = [];
        WM_agcc = [];
        GM_ctrl = [];
        GM_agcc = [];
        labels = label(idx_cap,NC);
        for i = 1:nframes_cap
            %Compute average with AgCC or controls only
            if strcmp(labels{i},'Ctrl')
                WM_ctrl = vertcat(WM_ctrl, avg_WM_truncated(i,:));
                GM_ctrl = vertcat(GM_ctrl, avg_GM_truncated(i,:));
            else
                WM_agcc = vertcat(WM_agcc, avg_WM_truncated(i,:));
                GM_agcc = vertcat(GM_agcc, avg_GM_truncated(i,:));
            end
        end
        WM_ctrl_list{l} = WM_ctrl;
        WM_agcc_list{l} = WM_agcc;
        GM_ctrl_list{l} = GM_ctrl;
        GM_agcc_list{l} = GM_agcc;
    
        avg_WM_both(l,:) = mean(avg_WM_truncated,1);
        avg_WM_ctrl(l,:) = mean(WM_ctrl,1);
        avg_WM_agcc(l,:) = mean(WM_agcc,1);
        avg_GM_both(l,:) = mean(avg_GM_truncated,1);
        avg_GM_ctrl(l,:) = mean(GM_ctrl,1);
        avg_GM_agcc(l,:) = mean(GM_agcc,1);
    end

    %% PCA 
    %For Inpainted signal 
    [XscoreIN,ncompIN] = apply_pca(active(idx_cap,:));
    %For WM signal 
    [XscoreWM,ncompWM] = apply_pca(active_maskedWM(idx_cap,:));
    %For GM signal 
    [XscoreGM,ncompGM] = apply_pca(active_maskedGM(idx_cap,:));

    %% SVM 100-fold
    fold = 100;
    accuracy_SVM_IN = zeros(1,fold);
    accuracy_SVM_WM = zeros(1,fold);
    accuracy_SVM_GM = zeros(1,fold);
    accuracy_SVM_WMpca = zeros(1,fold);
    accuracy_SVM_GMpca = zeros(1,fold);
    
    confusion_mat_IN = zeros(fold,2,2);
    confusion_mat_WM = zeros(fold,2,2);
    confusion_mat_GM = zeros(fold,2,2);
    confusion_mat_WMpca = zeros(fold,2,2);
    confusion_mat_GMpca = zeros(fold,2,2);
    
    %% Optimize hyperparameters
    
%     if ~exist(fullfile(savePath,strcat('SVMModel_linear',int2str(l),'.mat')),'file')
%         % Optimize the hyperparameters on whole dataset
%         SVMModel_IN = fitcsvm(XscoreIN(:,1:ncompIN),labels,...
%             'Standardize',true,'KernelFunction',char('linear'));%,...
%             %'OptimizeHyperparameters',{'BoxConstraint'});%,'KernelScale','KernelFunction'});
%         SVMModel_WMpca = fitcsvm(XscoreWM(:,1:ncompWM),labels,...
%             'Standardize',true,'KernelFunction',char('linear'));%,...
%             %'OptimizeHyperparameters',{'BoxConstraint'});%,'KernelScale','KernelFunction'});
%         SVMModel_GMpca = fitcsvm(XscoreGM(:,1:ncompGM),labels,...
%             'Standardize',true,'KernelFunction',char('linear'));%,...
%             %'OptimizeHyperparameters',{'BoxConstraint'});%,'KernelScale','KernelFunction'});
% 
%         save(fullfile(savePath,strcat('SVMModel_linear',int2str(l),'.mat')),'SVMModel_IN',...
%             'SVMModel_WMpca', 'SVMModel_GMpca');
%     else 
%         SVMmodels = load(fullfile(savePath,strcat('SVMModel_linear',int2str(l),'.mat')));
%         SVMModel_IN = SVMmodels.SVMModel_IN;
%         SVMModel_WMpca = SVMmodels.SVMModel_WMpca;
%         SVMModel_GMpca = SVMmodels.SVMModel_GMpca;
%     end
    
%     [Cin,gIN,kerIN] = get_hyperparameters(SVMModel_IN);
%     [Cwm,gWM,kerWM] = get_hyperparameters(SVMModel_WMpca);
%     [Cgm,gGM,kerGM] = get_hyperparameters(SVMModel_GMpca);

   %% Default Hypeparameters
    Cin = 1;
    gIN = 1;
    Cwm = 1;
    gWM = 1;
    Cgm = 1;
    gGM = 1;

    for i =1:fold
        %Build train and test sets
        disp(strcat('SVM fold ',int2str(i)))
        [idx_train,idx_test, label_train, label_test] = ...
        build_train_test_set(idx_cap,subjects_cap,i,NC,numsel_vec);

        %Build SVM
        SVMModel_IN = fitcsvm(XscoreIN(idx_train,1:ncompIN),label_train,...
            'Standardize',true,'KernelFunction',char('linear'),...
            'BoxConstraint',Cin,'KernelScale',gIN);
        SVMModel_WM = fitcsvm(avg_WM_truncated(idx_train,:),label_train,'Standardize',true);
        SVMModel_GM = fitcsvm(avg_GM_truncated(idx_train,:),label_train,'Standardize',true);
        SVMModel_WMpca = fitcsvm(XscoreWM(idx_train,1:ncompWM),label_train,...
            'Standardize',true,'KernelFunction',char('linear'),...
            'BoxConstraint',Cwm,'KernelScale',gWM);
        SVMModel_GMpca = fitcsvm(XscoreGM(idx_train,1:ncompGM),label_train,...
            'Standardize',true,'KernelFunction',char('linear'),...
            'BoxConstraint',Cgm,'KernelScale',gGM);

        %Evaluate SVM
        predictlabel_IN = predict(SVMModel_IN,XscoreIN(idx_test,1:ncompIN));
        accuracy_SVM_IN(i) = evaluate_SVM(predictlabel_IN,label_test);
        predictlabel_WM = predict(SVMModel_WM,avg_WM_truncated(idx_test,:));
        accuracy_SVM_WM(i) = evaluate_SVM(predictlabel_WM,label_test);
        predictlabel_GM = predict(SVMModel_GM,avg_GM_truncated(idx_test,:));
        accuracy_SVM_GM(i) = evaluate_SVM(predictlabel_GM,label_test); 
        predictlabel_WMpca = predict(SVMModel_WMpca,XscoreWM(idx_test,1:ncompWM));
        accuracy_SVM_WMpca(i) = evaluate_SVM(predictlabel_WMpca,label_test); 
        predictlabel_GMpca = predict(SVMModel_GMpca,XscoreGM(idx_test,1:ncompGM));
        accuracy_SVM_GMpca(i) = evaluate_SVM(predictlabel_GMpca,label_test); 
        
        confusion_mat_IN(i,:,:) = confusionmat(label_test,predictlabel_IN,'Order',{'AgCC','Ctrl'});
        confusion_mat_WM(i,:,:)= confusionmat(label_test,predictlabel_WM,'Order',{'AgCC','Ctrl'});
        confusion_mat_GM(i,:,:) = confusionmat(label_test,predictlabel_GM,'Order',{'AgCC','Ctrl'});
        confusion_mat_WMpca(i,:,:) = confusionmat(label_test,predictlabel_WMpca,'Order',{'AgCC','Ctrl'});
        confusion_mat_GMpca(i,:,:) = confusionmat(label_test,predictlabel_GMpca,'Order',{'AgCC','Ctrl'});

%         save(fullfile(savePath,strcat('predictlabel_cap_',int2str(l),...
%             '_fold',int2str(i),'.mat')),'predictlabel_IN',...
%             'predictlabel_WM','predictlabel_AH', 'predictlabel_WMpca',...
%             'predictlabel_GMpca','label_test'); 
    end
%     save(fullfile(savePath,strcat('accuracies_cap',int2str(l),'.mat')),...
%         'accuracy_SVM_AH', 'accuracy_SVM_WM');

    confusion_mat_IN_sum(l,:,:) = sum(confusion_mat_IN,1);
    confusion_mat_WM_sum(l,:,:) = sum(confusion_mat_WM,1);
    confusion_mat_GM_sum(l,:,:) = sum(confusion_mat_GM,1);
    confusion_mat_WMpca_sum(l,:,:) = sum(confusion_mat_WMpca,1);
    confusion_mat_GMpca_sum(l,:,:) = sum(confusion_mat_GMpca,1);

    mean_acc(1,l) = mean(accuracy_SVM_IN)*100;
    mean_acc(2,l) = mean(accuracy_SVM_WM)*100; 
    mean_acc(3,l) = mean(accuracy_SVM_GM)*100;
    mean_acc(5,l) = mean(accuracy_SVM_WMpca)*100;
    mean_acc(6,l) = mean(accuracy_SVM_GMpca)*100;
end 
save(fullfile(savePath,'SVM_class_results.mat'),'mean_acc');
save(fullfile(savePath,'average_per_cap.mat'),'avg_WM_both','avg_WM_agcc',...
    'avg_WM_ctrl','avg_GM_both','avg_GM_ctrl','avg_GM_agcc','WM_ctrl_list',...
    'WM_agcc_list','GM_ctrl_list','GM_agcc_list');
% save(fullfile(savePath,'average_in_cap.mat'),'avg_capWM','avg_capYEO');
save(fullfile(savePath,'confusion_matrices.mat'), 'confusion_mat_IN_sum',...
    'confusion_mat_WM_sum', 'confusion_mat_GM_sum', 'confusion_mat_WMpca_sum',...
    'confusion_mat_GMpca_sum');

function frame_out = apply_mask(mask,indices,frame_in)
    % Return the frame signal within the mask only
    % frame, mask are 1D vectors
    % indices is a 1D vector comprising the indices that we consider active to compute CAP
    mask = mask(indices);
    frame_out = frame_in .* mask';
end

function [Xscore,ncomp]= apply_pca(X)
%For GM signal
    [~,Xscore,~,~,explained,~] = pca(X);
    %nbr of components to take to explain 75% of the variance
    ncomp=1;
    while (sum(explained(1:ncomp))<75)
      ncomp= ncomp+1;  
    end
end

function [C,g,ker] = get_hyperparameters(model)
C = model.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint;
%g = model.HyperparameterOptimizationResults.XAtMinObjective.KernelScale;
g=1;
if isnan(g)
    g = 1; %default value of KernelScale
end
ker = 0;
%ker = char(model.HyperparameterOptimizationResults.XAtMinObjective.KernelFunction);
end

