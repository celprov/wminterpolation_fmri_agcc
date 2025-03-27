function subjects = dico_subject_category(listPath, dataPath)
%--------------------------------------------------------------------------
% Created by : Celine Provins (06.2020)
%
% Convert excel sheet into a cell of structure that pairs the subject names
% with their category
%
% INPUT
% listPath : foldername+filename of the excel sheet that repertoriates the
%   subject names with their category
% dataPath : path towards the folder containing the subject folders
% 
% OUTPUT
% subjects : cell of structure. One structure correspond to one subject and
%   contains the pairing between the name of the participant and its
%   category (AgCC partial or complete)
%--------------------------------------------------------------------------
    list_AgCC = readtable(listPath);
    list_name = list_AgCC.Var1;
    list_cat = list_AgCC.Var2;
    dirs = dir(fullfile(dataPath, 's*'));
    subjects = cell(size(dirs));
    %structure pairs the name of the participant present in the dataset
    %with their category (AgCC partial or complete)
    for i = 1:size(dirs,1)
        subject = struct();
        subject.name = dirs(i).name;
        for j = 1:length(list_name)
            if strcmp(list_name{j},dirs(i).name)
                subject.cat = list_cat{j};
            end
        end
        subjects{i} = subject; 
    end
end