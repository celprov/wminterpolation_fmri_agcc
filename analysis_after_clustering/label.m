function label_out = label(indices,NC)
%--------------------------------------------------------------------------
% Created by : Celine Provins (05.2020)
%
% Find out the label of each frame based on its index
%
% INPUT
% indices is a vector containing the index of the frames
% NC is the number of frames from control subjects
%
% OUTPUT 
% label_out : cell containing the label of each frame ('Ctrl' or 'AgCC')
%-------------------------------------------------------------------------- 
    
    label_out = cell(size(indices));
    for i = 1:length(indices)
        if indices(i)<=NC
            label_out{i}='Ctrl';
        else
            label_out{i}='AgCC';
        end       
    end
end  