function CAP = getCAPs(active_data,k,k_idx)
%--------------------------------------------------------------------------
% Created by : Zhiwei Huang (10.2019)
% Modified by : Celine Provins (04.2020)
%
% Build the CAPs by averaging the frames within clusters
% 
% INPUT
% active_data : matrix containing active frames
% k : integer corresponding to the number of clusters
% k_idx : vector containing the CAP membership number of each active frame
%
% OUTPUT
% CAP : matrix containing the built CAPs
%--------------------------------------------------------------------------
    for l = 1:k
        idx_l = find(k_idx==l);
        CAP(l,:) = mean(active_data(idx_l,:),1);
    end
end