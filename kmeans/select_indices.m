function [indices] = select_indices(active_data,threshold)
%--------------------------------------------------------------------------
% Created by : Zhiwei Huang (10.2019)
% 
% Select the voxels that are non-zero in at least the percentage of frames
% determined by the threshold
%
% INPUT
% active_data : matrix containing the active frames
% threshold : float [0,1] 
%
% OUTPUT
% indices : vector containing the linear coordinates of the voxel non-zero
%   in a certain percentage of the frames
%--------------------------------------------------------------------------
    active_data(find(active_data~=0))=1;
    projection = sum(active_data,1)/size(active_data,1);
    indices = find(projection>=threshold);
end