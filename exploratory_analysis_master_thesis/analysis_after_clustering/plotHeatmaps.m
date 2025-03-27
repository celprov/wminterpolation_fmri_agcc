% Created by : Celine Provins (06.2020)
% Uses the files obtained and saved by 'SVM_feature_2.m' and
% 'SVM_feature_1.m'
close all; clear all;

WholeBrain = 1; %work with WholeBrain CAPs (1) or PCC-seed CAPs (0)

%% Plot the heatmaps
celinePath = 'J:\Anjali_Diffusion_Pipeline\Celine';  
% celinePath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
codePath = fullfile(celinePath,'analysis_after_clustering');
savePathRFE = fullfile(codePath,'ResultsRFE_withACPC');
savePathClass = fullfile(codePath, 'ResultsClassification_withACPC');
statPath = fullfile(codePath,'CAPdistribution');
savePath = fullfile(codePath, 'Heatmaps','withACPC');
if WholeBrain
    statPath = fullfile(statPath,'WholeBrain');
    savePathRFE = fullfile(savePathRFE,'WholeBrain');
    savePathClass = fullfile(savePathClass,'WholeBrain');
    savePath = fullfile(savePath, 'WholeBrain');  
end

% load(fullfile(statPath,'statistics_avg'));
% load(fullfile(savePathClass,'average_per_cap.mat'));
% load(fullfile(savePathClass,'average_in_cap.mat'));
load(fullfile(savePathRFE,'RFE_results_default.mat'));
load(fullfile(savePathRFE,'accuracy_opt_default'));
load(fullfile(savePathRFE,'region_nameWM.mat'));

[signifWM, signifWM_top4, signifWM_top3] = build_ranking(ftRankWM);
[signifGM, signifGM_top4, signifGM_top3] = build_ranking(ftRankGM);

region_nameYEO =  {'VisCent', 'VisPeriph', 'SomotA', 'SomotB',...
                 'DorsAttA','DorsAttB',...
                 'VenAttn', 'Sal', 'LimbTemp', 'LimbOFC',...
                 'CtrlC', 'CtrlA', 'CtrlB', 'Temp/Par', 'DefC',...
                 'DefA', 'DefB'};
             
region_nameAH =  {'l-amPFC','l-PCC','dMPFC','l-TPJ','l-LTC','l-TempP','vMPFC',...
    'l-pIPL','l-Rsp','l-PHC','l-HF','r-amPFC','r-PCC','r-TPJ','r-LTC',...
    'r-TempP','r-pIPL','r-Rsp','r-PHC','r-HF' };

if WholeBrain
    % exclude CAP2,5 from analysis for WM
    signifWM([2,5],:)= [];
    signifWM_top4([2,5],:)= [];
    
     %exclude CAP5 from analysis for GM
    signifGM(5,:)= [];
    signifGM_top4(5,:)= [];
else
    % exclude CAP4,5,7 from analysis for WM
    signifWM([4,5,7],:)= [];
    signifWM_top4([4,5,7],:)= [];

    %exclude CAP1,4,9 and all from analysis for GM
    signifGM([1,4,9],:)= [];
    signifGM_top4([1,4,9],:)= [];
end
%% WM plots
% 
% figure 
% % subplot(2,2,1)
% yvalues = {'CAP-1','CAP-2','CAP-3','CAP-4','CAP-5','CAP-6','CAP-7','CAP-8'};
% h = heatmap(region_nameWM,yvalues,avg_WM_both,'CellLabelFormat','%07.4f');
% h.Colormap = jet;
% h.FontSize = 18;
% title('Average WM signal intensity')
% xlabel('WM bundles')
% ylabel('Number of CAPs')
% saveas(gcf,fullfile(savePath,strcat('heatmap_avgWM.png')));
% 
if WholeBrain
    yvalues = {'CAP-1','CAP-3','CAP-4','CAP-6','CAP-7','CAP-8','All'};
else
    yvalues = {'CAP-1','CAP-2','CAP-3','CAP-6','CAP-8','All'};
end
figure
subplot(2,2,1)

h = heatmap(region_nameWM,yvalues,signifWM,'CellLabelFormat','%07.4f');
h.Colormap = flipud(hot);
h.FontSize = 13;
title('A. WM feature ranking')
xlabel('WM bundles')
ylabel('Number of CAPs')

subplot(2,2,2)
h = heatmap(region_nameWM,yvalues,signifWM_top4,'CellLabelFormat','%07.4f');
h.Colormap = flipud(hot);
h.GridVisible = 'off';
h.FontSize = 13;
title('B. WM feature ranking - Top 4')
xlabel('WM bundles')
ylabel('Number of CAPs')

% subplot(2,2,3)
% yvalues ={'CAP-1','CAP-2','CAP-3','CAP-4','CAP-5','CAP-6','CAP-7','CAP-8'};
% h = heatmap(region_nameWM,yvalues,hWM,'CellLabelFormat','%07.4f');
% % h.Colormap = jet;
% h.FontSize = 13;
% title('WM significance')
% xlabel('WM bundles')
% ylabel('Number of CAPs')
% saveas(gcf,fullfile(savePath,'heatmapsWMsignificance.png'))

% subplot(2,2,3)
% yvalues = {'CAP-1','CAP-2','CAP-3','CAP-4','CAP-5','CAP-6','CAP-7','CAP-8'};
% h = heatmap(region_nameWM,yvalues,avg_WM_ctrl,'CellLabelFormat','%07.4f');
% h.Colormap = jet;
% h.FontSize = 13;
% title('C. Average WM signal intensity - Control')
% xlabel('WM bundles')
% ylabel('Number of CAPs')
% 
% subplot(2,2,4)
% yvalues = {'CAP-1','CAP-2','CAP-3','CAP-4','CAP-5','CAP-6','CAP-7','CAP-8'};
% h = heatmap(region_nameWM,yvalues,avg_WM_agcc,'CellLabelFormat','%07.4f');
% h.Colormap = jet;
% h.ColorLimits = [min(avg_WM_ctrl(:)),max(avg_WM_ctrl(:))];
% h.FontSize = 13;
% title('D. Average WM signal intensity - AgCC')
% xlabel('WM bundles')
% ylabel('Number of CAPs')
% % saveas(gcf,fullfile(savePath,strcat('heatmapsWM.png')));

%% GM plots
% figure 
% % subplot(2,2,1)
% yvalues = {'CAP-1','CAP-2','CAP-3','CAP-4','CAP-5','CAP-6','CAP-7','CAP-8'};
% h = heatmap(region_nameAH,yvalues,avg_AH_both,'CellLabelFormat','%07.4f');
% h.Colormap = jet;
% h.FontSize = 18;
% title('Average GM signal intensity')
% xlabel('GM regions')
% ylabel('Number of CAPs')
% % saveas(gcf,fullfile(savePath,strcat('heatmap_avgAH.png')));
if WholeBrain
    yvalues = {'CAP-1','CAP-2','CAP-3','CAP-4','CAP-6','CAP-7','CAP-8','All'};
    region_nameGM = region_nameYEO;
else
    yvalues = {'CAP-2','CAP-3','CAP-5','CAP-6','CAP-7','CAP-8'};
    region_nameGM = region_nameAH;
end 
figure 
subplot(2,2,1)
h = heatmap(region_nameGM,yvalues,signifGM,'CellLabelFormat','%07.4f');
h.Colormap = flipud(hot);
h.GridVisible = 'off';
h.FontSize = 13;
title('A. GM feature ranking')
xlabel('GM regions')
ylabel('Number of CAPs')

subplot(2,2,2)
h = heatmap(region_nameGM,yvalues,signifGM_top4,'CellLabelFormat','%07.4f');
h.Colormap = flipud(hot);
h.GridVisible = 'off';
h.FontSize = 13;
title('B. GM feature ranking - Top 4')
xlabel('GM regions')
ylabel('Number of CAPs')

subplot(2,2,3)
yvalues ={'CAP-1','CAP-2','CAP-3','CAP-4','CAP-5','CAP-6','CAP-7','CAP-8'};
h = heatmap(region_nameGM,yvalues,hGM,'CellLabelFormat','%07.4f');
% h.Colormap = jet;
h.FontSize = 13;
title('GM significance')
xlabel('GM regions')
ylabel('Number of CAPs')
% saveas(gcf,fullfile(savePath,'heatmapsGMsignificance.png'))
 
% subplot(2,2,3)
% yvalues = {'CAP-1','CAP-2','CAP-3','CAP-4','CAP-5','CAP-6','CAP-7','CAP-8'};
% h = heatmap(region_nameGM,yvalues,avg_GM_ctrl,'CellLabelFormat','%07.4f');
% h.Colormap = jet;
% h.FontSize = 13;
% h.ColorLimits = [min(avg_GM_agcc(:)),max(avg_GM_agcc(:))];
% title('C. Average GM signal intensity - Control')
% xlabel('GM regions')
% ylabel('Number of CAPs')
% 
% subplot(2,2,4)
% yvalues = {'CAP-1','CAP-2','CAP-3','CAP-4','CAP-5','CAP-6','CAP-7','CAP-8'};
% h = heatmap(region_nameGM,yvalues,avg_GM_agcc,'CellLabelFormat','%07.4f');
% h.Colormap = jet;
% h.FontSize = 13;
% title('D. Average GM signal intensity - AgCC')
% xlabel('GM regions')
% ylabel('Number of CAPs')
% saveas(gcf,fullfile(savePath,strcat('heatmapsGM.png')));

%% Plot the accuracy with respect to threshold
figure
subplot(1,2,2)
thresh = [1,0.9,0.8,0.7,0.6,0.5];
for l = 1:9
    p = plot(thresh,mean_acc_GM(:,l),'-o','LineWidth',1.5);
    if l == 9
        p.Color = [0.04,0.30,0.08];
    elseif l==8
        p.Color = [1.00,0.07,0.65];
    end
    hold on 
end
legend('CAP1','CAP2','CAP3','CAP4','CAP5','CAP6','CAP7','CAP8','All','Location','northeastoutside')
xlabel('Correlation threshold')
ylabel('Accuracy (%)')
set(gca,'Fontsize',16)
xlim([0.5,1])
title(sprintf('B. Classification accuracy of SVM\nwith GM signal as features'))

subplot(1,2,1)
for l = 1:9
    p = plot(thresh,mean_acc_WM(:,l),'-o','LineWidth',1.5);
    if l == 9
        p.Color = [0.04,0.30,0.08];
    elseif l == 8
        p.Color = [1.00,0.07,0.65];
    end
    hold on 
end
legend('CAP1','CAP2','CAP3','CAP4','CAP5','CAP6','CAP7','CAP8','All','Location','northeastoutside')
xlabel('Correlation threshold')
ylabel('Accuracy (%)')
set(gca,'Fontsize',16)
xlim([0.5,1])
title(sprintf('A. Classification accuracy of SVM\nwith WM signal as features'))

function [signif,signif_top4, signif_top3] = build_ranking(ftRank)
    % Build the ranking 
    % ftRank is a vector containing the ranked features
    k = size(ftRank,1);
    nbr_region = size(ftRank,2);
    
    signif = zeros(k,nbr_region);
    for l = 1:k 
        rank = nbr_region;
        for i = 1:nbr_region
            signif(l,ftRank(l,i)) = rank;
            rank = rank-1;
        end
    end
    
    % Build the ranking for the top 4 features
    signif_top4 = zeros(k,nbr_region);
    for l = 1:k
        idx4 = find(signif(l,:) == nbr_region);
        idx3 = find(signif(l,:) == nbr_region-1);
        idx2 = find(signif(l,:) == nbr_region-2);
        idx1 = find(signif(l,:) == nbr_region-3);
        signif_top4(l,idx4) = 4;
        signif_top4(l,idx3) = 3;
        signif_top4(l,idx2) = 2;
        signif_top4(l,idx1) = 1;
    end
    
     % Build the ranking for the top 3 features
    signif_top3 = zeros(k,nbr_region);
    for l = 1:k
        idx3 = find(signif(l,:) == nbr_region);
        idx2 = find(signif(l,:) == nbr_region-1);
        idx1 = find(signif(l,:) == nbr_region-2);
        signif_top3(l,idx3) = 3;
        signif_top3(l,idx2) = 2;
        signif_top3(l,idx1) = 1;
    end
end