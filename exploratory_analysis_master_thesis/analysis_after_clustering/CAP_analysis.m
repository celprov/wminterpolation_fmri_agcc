% Created by : Celine Provins (23.04.2020)
clear all; close all;
warning('off')

WholeBrain = 0; %work with WholeBrain CAPs (1) or PCC-seed CAPs (0)

%% Set up paths
% celinePath = 'J:\Anjali_Diffusion_Pipeline\Celine';  
celinePath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
codePath = fullfile(celinePath, 'analysis_after_clustering');
dataPath = fullfile(celinePath,'data');
agccPath = fullfile(dataPath,'RestingState');
ctrlPath = fullfile(dataPath, 'ControlsRS');
if WholeBrain 
    capPathBoth = fullfile(dataPath,'WholeBrainCAP_RS');
    capPathA = fullfile(agccPath,'WholeBrainCAP');
    capPathC = fullfile(ctrlPath,'WholeBrainCAP');
    savePath = fullfile(codePath,'CAPdistribution','WholeBrain');
    savePathClass = fullfile(codePath, 'ResultsClassification_withACPC','WholeBrain');
    savePathRFE = fullfile(codePath, 'ResultsRFE_withACPC','WholeBrain');
else
    capPathBoth = fullfile(dataPath,'CAP_RS_zscore');
    capPathA = fullfile(agccPath,'CAP');
    capPathC = fullfile(ctrlPath,'CAP');
    savePath = fullfile(codePath,'CAPdistribution');
    savePathClass = fullfile(codePath, 'ResultsClassification_withACPC');
    savePathRFE = fullfile(codePath, 'ResultsRFE_withACPC');
end
listPath = fullfile(dataPath, 'list_partial_full_AgCC.csv');
maskPath = fullfile(celinePath,'Preprocessing-Hub','NuisanceRegression');
spmPath = fullfile(celinePath,'spm12');

addpath(genpath(spmPath));
addpath(genpath(codePath));

%% Load variables
load(fullfile(savePathClass,'average_per_cap.mat'));
load(fullfile(savePathRFE,'region_nameWM.mat'));

frame_check = load(fullfile(capPathBoth,'frame_check.mat'));
frame_GSR = frame_check.frame_GSR;

activeA = load(fullfile(capPathA,'active.mat'));
activeA_truncated=activeA.active_GSR_truncated;
activeC = load(fullfile(capPathC,'active.mat'));
activeC_truncated=activeC.active_GSR_truncated;

% cap = load(fullfile(capPath,'CAP.mat'));
% cap_gsr = cap.CAP_both_GSR;
% idx = cap.idx;

%% Define constants

%Numbers for Working Memory
% numsel = 57; %nbr of frame selected
% nc = 16; %nbr control
% na = 7; %nbr agcc
% nap = 2;%nbr agcc partial
% nac = 5;%nbr agcc complete
 
%Numbers for Resting State
numsel = 30; %nbr of frame selected
nc = 27; %nbr control
na = 15; %nbr agcc subjects
nap = 5;%nbr agcc partial subjects
nac = 10;%nbr agcc complete subjects
n = nc + nac + nap; %total nbr of subjects

NC = size(activeC_truncated,1); %nbr of frames from controls
NA = size(activeA_truncated,1); %nbr of frames from AgCC
nframes = NA + NC;

%the list of indices specific to each cap are separated by zeroes
delimiter_GSR = find(frame_GSR==0);
k = length(delimiter_GSR); %nbr of clusters

%Build the vector that indicates to which subject belong the frame i
numsel_vec = build_numsel_vec(WholeBrain,activeA,activeC,NC,nframes,numsel);

%Structure pairing subject name with partial or complete AgCC
AgCCsubjects = dico_subject_category(listPath, agccPath);
NAC = 0; %nbr of frames from complete AgCC
NAP = 0; %nbr of frames from partial AgCC
for i = NC+1:nframes
    m = numsel_vec(i)-nc;
    subject = AgCCsubjects{m};
    if strcmp(subject.cat, '1')
        NAC = NAC+1;
    end
    if strcmp(subject.cat, '2')
        NAP = NAP+1;
    end
end

%% Extracting the frames indices contributing to each CAP
list_idx_cap_nonordered = {1,k}; %necessary because CAP_RS and CAP_RS_zscore not same #CAP
list_idx_cap = {1,k};
for l = 1:k
    if l~=k
        list_idx_cap_nonordered{l} = frame_GSR(delimiter_GSR(l)+1:delimiter_GSR(l+1)-1);
    else
        list_idx_cap_nonordered{l} = frame_GSR(delimiter_GSR(l)+1:end);
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

%% Distribution of partial and complete AgCC as well as controls
% % the index contributing to the cap_i is identified from the partial or
% % complete category by linking its dividende to the name of the
% % corresponding subject
% % the indices <= 30 are coming from the first healthy subject
% % the indices 30 <= 60 are coming from the second healthy subject...
caps = {1,k};
for i= 1:k
    cap = struct();
    idx_cap = list_idx_cap{i};
    %Determine to which group each frame belongs
    idx_cap_partial_frame = [];
    idx_cap_complete_frame = [];
    idx_cap_ctrl_frame = [];
    idx_cap_partial_subject = [];
    idx_cap_complete_subject = [];
    idx_cap_ctrl_subject = [];
    for j = 1:length(idx_cap)
        idx = idx_cap(j);
        m = numsel_vec(idx);
        if m>nc
            subject = AgCCsubjects{m-nc};
            if strcmp(subject.cat, '1')
                idx_cap_complete_frame = [idx_cap_complete_frame idx];
                idx_cap_complete_subject = [idx_cap_complete_subject m];
            end
            if strcmp(subject.cat, '2')
                idx_cap_partial_frame = [idx_cap_partial_frame idx];
                idx_cap_partial_subject = [idx_cap_partial_subject m];
            end
        else
            idx_cap_ctrl_frame = [idx_cap_ctrl_frame idx]; 
            idx_cap_ctrl_subject = [idx_cap_ctrl_subject m]; 
        end
    end
    cap.partial_frame =  idx_cap_partial_frame;
    cap.complete_frame = idx_cap_complete_frame;
    cap.ctrl_frame = idx_cap_ctrl_frame;
    
    cap.partial_subject =  unique(idx_cap_partial_subject);
    cap.complete_subject = unique(idx_cap_complete_subject);
    cap.ctrl_subject = unique(idx_cap_ctrl_subject);
    
    caps{i}=cap;
end

dist_bar = zeros(k,2);
dist_bar_agcc = zeros(k,2);
data_bar = zeros(k,3);
for i = 1:k
    NAP_cap = length(caps{i}.partial_frame); %nbr of partial AgCC in CAP-i
    NAC_cap = length(caps{i}.complete_frame); %nbr of complete AgCC in CAP-i
    NC_cap = length(caps{i}.ctrl_frame); %nbr of control in CAP-i
    data_bar(i,1) = NC_cap;
    data_bar(i,2) = NAC_cap;
    data_bar(i,3) = NAP_cap;
    
    %weight to compensate difference in population number
    NC_cap = NC_cap / NC;
    NAC_cap = NAC_cap /NAC;
    NAP_cap = NAP_cap /NAP;
    
    dist_bar(i,1) = NC_cap/(NC_cap+NAC_cap+NAP_cap)*100;
    dist_bar(i,2) = NAC_cap/(NC_cap+NAC_cap+NAP_cap)*100;
    dist_bar(i,3) = NAP_cap/(NC_cap+NAC_cap+NAP_cap)*100;
end

data_bar(:,1) = data_bar(:,1)./NC.*100;
data_bar(:,2) = data_bar(:,2)./NAC.*100;
data_bar(:,3) = data_bar(:,3)./NAP.*100;

save(fullfile(savePath,'data_to_plot.mat'),'dist_bar','data_bar');

%% Plot CAP distribution - bar plot
% figure
% label_bar = categorical({'CAP1','CAP2','CAP3','CAP4','CAP5','CAP6',...
%     'CAP7','CAP8','All'});
% bar(label_bar,data_bar,'stacked')
% ylabel('Nbr of subjects')
% set(gca,'fontsize',13)
% legend('Control','complete AgCC', 'partial AgCC');

figure
subplot(1,4,1)
label_bar = categorical({'CAP1','CAP2','CAP3','CAP4','CAP5','CAP6',...
    'CAP7','CAP8'});
bar(label_bar,dist_bar,'stacked')
yline(33,'--','LineWidth',1.2,'Color','k')
yline(66,'--','LineWidth',1.2,'Color','k')
yticks([0,33,66,100])
ylim([0,100])
ylabel('CAPs percentage')
set(gca,'fontsize',17)
legend('Control','Complete AgCC', 'Partial AgCC','Location','south');
title('A. Distribution within CAP')

subplot(1,4,2)
label_bar = categorical({'CAP1','CAP2','CAP3','CAP4','CAP5','CAP6',...
    'CAP7','CAP8'});
bar(label_bar,data_bar(:,1),'FaceColor',[0 0.4470 0.7410])
ylabel('Nbr of frames contributing(%)')
ylim([0,35])
set(gca,'fontsize',17)
title('B. Controls');

subplot(1,4,3)
label_bar = categorical({'CAP1','CAP2','CAP3','CAP4','CAP5','CAP6',...
    'CAP7','CAP8'});
bar(label_bar,data_bar(:,2),'FaceColor',[0.8500 0.3250 0.0980])
ylabel('Nbr of frames contributing(%)')
ylim([0,35])
set(gca,'fontsize',17)
title('C. AgCC complete');

subplot(1,4,4)
label_bar = categorical({'CAP1','CAP2','CAP3','CAP4','CAP5','CAP6',...
    'CAP7','CAP8'});
bar(label_bar,data_bar(:,3),'FaceColor',[0.9290 0.6940 0.1250])
ylabel('Nbr of frames contributing (%)')
ylim([0,35])
set(gca,'fontsize',17)
title('D. AgCC partial');
saveas(gcf,fullfile(savePath,'CAPdist.png'));

% % subplot(1,3,3)
% % label_bar = categorical({'CAP1','CAP2','CAP3','CAP4','CAP5','CAP6',...
% %     'CAP7','CAP8'});
% % bar(label_bar,agcc_data,'FaceColor',[0.8500 0.3250 0.0980])
% % ylabel('Nbr of frames contributing (%)')
% % set(gca,'fontsize',16)
% % title('C. AgCC');
% 
% figure
% figure
% subplot(1,3,1)
% label_bar = categorical({'CAP1','CAP2','CAP3','CAP4','CAP5','CAP6',...
%     'CAP7','CAP8'});
% b = bar(label_bar,dist_bar_agcc,'stacked');
% b(1).FaceColor = [0.9290 0.6940 0.1250];
% b(2).FaceColor = [0.4940 0.1840 0.5560];
% yline(0.5,'--','LineWidth',1.2,'Color','k')
% yticks([0,0.25,0.5,0.75,1])
% ylabel('CAPs percentage')
% set(gca,'fontsize',16)
% legend('Complete AgCC','Partial AgCC');
% title('A. Distribution within CAP')
% 
% subplot(1,3,2)
% label_bar = categorical({'CAP1','CAP2','CAP3','CAP4','CAP5','CAP6',...
%     'CAP7','CAP8'});
% bar(label_bar,data_bar(:,2),'FaceColor',[0.9290 0.6940 0.1250])
% ylabel('Nbr of frames contributing(%)')
% ylim([0,40])
% set(gca,'fontsize',16)
% title('B. AgCC complete');
% 
% subplot(1,3,3)
% label_bar = categorical({'CAP1','CAP2','CAP3','CAP4','CAP5','CAP6',...
%     'CAP7','CAP8'});
% bar(label_bar,data_bar(:,3),'FaceColor',[0.4940 0.1840 0.5560])
% ylabel('Nbr of frames contributing (%)')
% ylim([0,40])
% set(gca,'fontsize',16)
% title('C. AgCC partial');
% 
% figure
% label_bar = categorical({'CAP1','CAP2','CAP3','CAP4','CAP5','CAP6',...
%     'CAP7','CAP8','All'});
% bar(label_bar,data_bar(:,2),'r')
% ylabel('Nbr of subjects')
% set(gca,'fontsize',13)
% legend('Complete AgCC');
% 
% figure
% label_bar = categorical({'CAP1','CAP2','CAP3','CAP4','CAP5','CAP6',...
%     'CAP7','CAP8','All'});
% bar(label_bar,data_bar(:,3),'y')
% ylabel('Nbr of subjects')
% set(gca,'fontsize',13)
% legend('Partial AgCC');
% saveas(gcf,fullfile(capPath,'distribution_bar.png'));

%% Investigate CAP occurence
% Compute CAP occurence
occurence=zeros(n,k);
for i = 1:n
    for l = 1:k
        %Extracting frames contributing to the CAP l
        idx_cap = list_idx_cap{l};
        %identify indices from the subject
        idx_sub = find(numsel_vec == i);
        % indices belonging to the subject and the CAP
        idx_sub_CAP = intersect(idx_sub,idx_cap);
%         idx_sub_CAP = find(idx_cap <= i*numsel & idx_cap>(i-1)*numsel);
        %compute CAP occurence in each subject
        occurence(i,l) = length(idx_sub_CAP)/length(idx_sub);
    end
end

% Investigate statistical difference between the occurences of each group
h_norm=zeros(2,k);
h=zeros(1,k);
p=zeros(1,k);
h3=zeros(1,k);
p3=zeros(1,k);
mean_ctrl = zeros(1,k);
mean_complete = zeros(1,k);
mean_partial = zeros(1,k);
stderror_ctrl = zeros(1,k);
stderror_complete = zeros(1,k);
stderror_partial = zeros(1,k);
for l = 1:k
    ctrl_group = caps{l}.ctrl_subject;
    partial_group = caps{l}.partial_subject;
    complete_group = caps{l}.complete_subject;
    mc = length(ctrl_group);
    mac = length(complete_group);
    map = length(partial_group);
    
    %Build the reference label vector 
    label = zeros(1,mc+mac+map);
    for i = 1:mc
    label(i) = 1; %control
    end
    for i = mc+1:mc+mac
    label(i) = 2; %complete AgCC
    end
    for i = mc+mac+1:mc+mac+map
    label(i) = 3; % partial AgCC
    end
    
    %check if occurrences are coming from a normal distribution
    h_norm(1,l) = adtest(occurence(1:nc,l));
    h_norm(2,l) = adtest(occurence(nc+1:n,l));
    %check if the samples come from same distribution with same median
    %Control vs AgCC - test with two group
    [p(l), h(l)] = ranksum(occurence(ctrl_group,l),...
        vertcat(occurence(complete_group,l), occurence(partial_group,l)));
    %Control vs complete AgCC vs partial AgCC - test with three groups
    p3(l) = kruskalwallis(...
        horzcat(occurence(ctrl_group,l)', occurence(complete_group,l)', occurence(partial_group,l)'),...
        label,'off');
    
    %Mean occurence
    mean_ctrl(l) = mean(occurence(ctrl_group,l),1)*100;
    mean_complete(l) = mean(occurence(complete_group,l),1)*100;
    mean_partial(l) = mean(occurence(partial_group,l),1)*100;
    %Standard error of the occurence
    stderror_ctrl(l) = std(occurence(ctrl_group,l),0,1) / sqrt(length(occurence(ctrl_group,l)))*100;
    stderror_complete(l) = std(occurence(complete_group,l),0,1) / sqrt(length(occurence(complete_group,l)))*100;
    stderror_partial(l) = std(occurence(partial_group,l),0,1) / sqrt(length(occurence(partial_group,l)))*100;
end

save(fullfile(savePath,'statistics_occurence.mat'),'h_norm','p','h','p3');

%Plot occurence results
cap = 1:8;
figure
errorbar(cap,mean_ctrl,stderror_ctrl,'o','LineWidth',1.2,'CapSize',10)
hold on 
errorbar(cap,mean_complete,stderror_complete,'^','LineWidth',1.2,'CapSize',10)
e = errorbar(cap,mean_partial,stderror_partial,'x','LineWidth',1.2,'CapSize',10);
e.Color = [0.4940 0.1840 0.5560];
xlim([0.5,8.5])
set(gca,'FontSize',16)
XTicks = [1,2,3,4,5,6,7,8];
% plot(2,2.7,'*k','LineWidth',1) %plot the asteriks indicating significance
% plot(3,1.55,'*k','LineWidth',1)
% ylim([0,3])
legend('Control','Complete AgCC', 'Partial AgCC','Location','north')
ylabel('Occurence (%)')
xlabel('Number of CAP')
saveas(gcf,fullfile(savePath,'CAPoccurence.png'));

%% Statistical tests on average within CAP
nbr_regionWM = size(WM_ctrl_list{1},2);
nbr_regionGM = size(GM_ctrl_list{1},2);
pWM = zeros(k,nbr_regionWM);
hWM = zeros(k,nbr_regionWM);
pGM = zeros(k,nbr_regionGM);
hGM = zeros(k,nbr_regionGM);
h_norm_avgWM = zeros(2,k,nbr_regionWM);
h_norm_avgGM = zeros(2,k,nbr_regionGM);
for l = 1:k
    WM_ctrl = WM_ctrl_list{l};
    WM_agcc = WM_agcc_list{l};
    GM_ctrl = GM_ctrl_list{l};
    GM_agcc = GM_agcc_list{l};
    
    for i = 1 : nbr_regionWM
        %check if occurrences are coming from a normal distribution
        h_norm_avgWM(1,l,i) = adtest(WM_ctrl(:,i));
        h_norm_avgWM(2,l,i) = adtest(WM_agcc(:,i));
        
        %check if the samples come from same distribution with median of
        %agcc higher than median of ctrl
        %Control vs AgCC - test with two group
        [pWM(l,i), hWM(l,i)] = ranksum(WM_ctrl(:,i),WM_agcc(:,i),'tail','left');
        
    end
    
    for i = 1 : nbr_regionGM
         h_norm_avgGM(1,l,i) = adtest(GM_ctrl(:,i));
         h_norm_avgGM(2,l,i) = adtest(GM_agcc(:,i));
         [pGM(l,i), hGM(l,i)] = ranksum(GM_ctrl(:,i),GM_agcc(:,i),'tail','left');
    end
end

save(fullfile(savePath,'statistics_avg.mat'),'h_norm_avgWM','h_norm_avgGM',...
    'pWM','hWM','pGM','hGM');

