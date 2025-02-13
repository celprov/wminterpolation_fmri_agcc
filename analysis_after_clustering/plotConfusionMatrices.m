%Created by : Celine Provins (06.2020)
% Uses the files obtained and saved by 'SVM_classification.m'
close all; clear all;

WholeBrain = 0; %work with WholeBrain CAPs (1) or PCC-seed CAPs (0)

%% Plot the confusion matrices
% celinePath = 'J:\Anjali_Diffusion_Pipeline\Celine';  
celinePath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
codePath = fullfile(celinePath,'analysis_after_clustering');
savePathClass = fullfile(codePath, 'ResultsClassification_withACPC');
savePath = fullfile(codePath, 'ConfusionMatrices','WithACPC');
if WholeBrain
    savePathClass = fullfile(savePathClass,'WholeBrain');
    savePath = fullfile(savePath,'WholeBrain');
end

load(fullfile(savePathClass,'confusion_matrices.mat'));

plot_confusion_mat(confusion_mat_IN_sum,'Inpainted',savePath);
plot_confusion_mat(confusion_mat_WM_sum, 'Average_WM',savePath);
plot_confusion_mat(confusion_mat_GM_sum, 'Average_GM',savePath);
plot_confusion_mat(confusion_mat_WMpca_sum, 'WM_mask',savePath);
plot_confusion_mat(confusion_mat_GMpca_sum, 'GM_mask',savePath);

function plot_confusion_mat(confusion_mat,titre,savePath)
    k = size(confusion_mat,1);
    figure
    title(titre)
    for l = 1:k
        C = confusionchart(squeeze(confusion_mat(l,:,:)),{'Ctrl','AgCC'},...
            'Normalization', 'row-normalized', 'DiagonalColor',...
            [0.6350 0.0780 0.1840], 'OffDiagonalColor',[0.6350 0.0780 0.1840]);
        C.FontSize = 50;
        C.XLabel ='';
        C.YLabel ='';
        saveas(gcf,fullfile(savePath,strcat('conf_mat_',titre,'CAP',int2str(l),'.png'))) 
    end  
end