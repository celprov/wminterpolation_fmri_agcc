function [k_opt] = find_optimal_k(sample_data,k_range,capPath)
%--------------------------------------------------------------------------
% Created by : Zhiwei Huang (10.2019)
%
% Find optimal number of clusters based on silhouette coefficient
%
% INPUT
% sample_data : matrix containing active frames
% k_range : integer corresponding to the range of the number of CAPs k
% capPath : string giving the path to save optimal k search results
%
% OUTPUT
% k_opt : the optimal number of clusters corresponding with k presenting
%   the highest silhouette coefficient
%--------------------------------------------------------------------------
    for k = 1:length(k_range)
        idx_k = kmeans(sample_data,k_range(k)); 
        s = silhouette(sample_data,idx_k,'Euclidean');
        s_mean(k) = mean(s);
    end
    close all
    plot(k_range,s_mean)
    hold on
    title('Mean Silhouette Value over K')
    xlabel('Number of Clusters')
    ylabel('Silhouette Value')
    saveas(gcf,fullfile(capPath,'k_optimal.png'))
    %optimal k corresponds to the k having maximal silhouette coefficient
    k_opt_idx = find(s_mean==max(s_mean));
    k_opt = k_range(k_opt_idx(1));
end