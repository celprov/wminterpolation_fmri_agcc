function [idx_train,idx_test, label_train, label_test] = ...
    build_train_test_set(idx_cap,subjects_cap,fold,NC,numsel_vec)
%--------------------------------------------------------------------------
% Created by : Celine Provins (05.2020)
%
% Separate the indices contained in idx_cap into a train and a test set.
% The frames corresponding to one particular subject (determined by fold)
% are removed from the train set and added in the test set.
% 
% INPUT
% idx_cap is a 1D vector containing the indices of the frames contributing
%   to the CAP
% subjects_cap is a 1D vector containing the index of the subject
%   contributing to the CAP
% fold is an int corresponding to fold number
% NC is the number of controls frames
% numsel is the number of frames selected per subject
%
% OUTPUT
% idx_train : indices of the frames put in the train set
% idx_test : indices of the frames put in the test set
% label_train : labels corresponding to the frames in the train set
% label_test : labels corresponding to the frames in the test set
%--------------------------------------------------------------------------
 
    %Use one subject at each fold as test set
    nsub=length(subjects_cap);
    nframes_cap = length(idx_cap);
    u = mod(fold,nsub);
    if u ==0 %subject 0 doesn't exist
        u = nsub;
    end
    s = subjects_cap(u); 
    % remove indices corresponding to that subject from training set
    % and add it in test set
    trainset = idx_cap;
    idx_train = 1:nframes_cap;
    testset = [];
    idx_test = [];
    for j = 1:nframes_cap
        if numsel_vec(idx_cap(j)) == s
            testset = [testset idx_cap(j)];
            idx_test = [idx_test j];
            trainset(j) = 0; 
            %set to zero and then remove because if remove straight
            %away problem of dimensionality
            idx_train(j) = 0;
        end
    end
    trainset = nonzeros(trainset);
    idx_train = nonzeros(idx_train);

%   save(fullfile(savePath,strcat('idx_set_cap',int2str(l),...
%             '_fold',int2str(fold),'.mat')),'idx_train','idx_test');
%   save(fullfile(savePath,strcat('train_test_set_cap',int2str(l),...
%             '_fold',int2str(fold),'.mat')),'trainset','testset');

    % Find the label corresponding to each frame
    label_train = label(trainset,NC);
    label_test = label(testset,NC);

%   save(fullfile(savePath,strcat('label_cap',int2str(l),...
%    '_fold',int2str(i),'.mat')),'label_test','label_train');
end