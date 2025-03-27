function numsel_vec = build_numsel_vec(WholeBrain,activeA,activeC,NC,nframes,numsel)
%--------------------------------------------------------------------------
% Created by : Celine Provins (06.2020)
%
% Build the vector that indicates to which subject belong the frames
% Necessary because the number of frame selected per patient can be
% different in the case of the WholeBrain
%
% INPUT
% WholeBrain : bool that precise whether we work with Whole Brain CAP (1)
%    or PCC-seed CAPs (0)
% activeA : struct where idx_active is stored for AgCC patients
% activeC : struct where idx_active is stored for control patients
% NC : number of control active frames 
% nframes : total number of frames in dataset
% numsel : number of active frames selected per subject
%
% OUTPUT
% numsel_vec : vector (1xnframes) that containes to which subject belongs
%   each frame
%--------------------------------------------------------------------------
if WholeBrain 
    idx_activeC = activeC.idx_active; 
    idx_activeA = activeA.idx_active;
    nc = length(idx_activeC); %nbr of control subjects
    na = length(idx_activeA); %nbr of AgCC subjects
    numsel_vec = zeros(1,nframes);
    numsel_sub_1 = 1;
    for j = 1:nc
       numsel_sub_2 = length(idx_activeC{j});
       numsel_vec(numsel_sub_1:numsel_sub_1+numsel_sub_2-1) = j;
       numsel_sub_1 = numsel_sub_1+numsel_sub_2;
    end
    numsel_sub_1 = 1;
    for j = 1:na
       numsel_sub_2 = length(idx_activeA{j});
       numsel_vec(numsel_sub_1+NC:numsel_sub_1+numsel_sub_2+NC-1) = j+nc;
       numsel_sub_1 = numsel_sub_1+numsel_sub_2;
    end    
else
    nc = NC/numsel;%nbr of control subjects
    na = (nframes-NC)/numsel; %nbr of AgCC subjects
    for j = 0:nc-1
        numsel_vec(1+numsel*j:numsel*(j+1)) = j+1;
    end
    for j = 0:na-1
        numsel_vec(1+numsel*j+NC:numsel*(j+1)+NC) = j+1+nc;
    end
end