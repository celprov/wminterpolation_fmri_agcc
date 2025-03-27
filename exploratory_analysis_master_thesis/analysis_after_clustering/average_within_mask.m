function avg = average_within_mask(list_masks,indices,frame)
%--------------------------------------------------------------------------
% Created by : Celine Provins (05.2020)
%
% Compute the average signal of the frame within the mask
%
% INPUT
% list_masks : cell (1 x #regions) containing the mask put in 1D of each region
% frame : vector containing the frame to average within the mask
% indices : 1D vector comprising the linear indices of the frames
%   considered active
%
% OUTPUT
% avg : vector (1 x #regions) containing the average of the frame within
%   the atlas
%--------------------------------------------------------------------------
    avg = zeros(1,length(list_masks));
    for j = 1:length(list_masks)
        mask = list_masks{j};
        indices_mask = indices;
        %compute the average of each frame within mask
        for l = 1:length(indices)
            %selecting if the index belongs to the mask or not
            indices_mask(l)=indices(l)*mask(indices(l));
        end
        avg(j) = mean(frame(find(indices_mask)));
    end  
end