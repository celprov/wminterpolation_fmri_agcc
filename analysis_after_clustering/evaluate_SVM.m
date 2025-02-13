function accuracy = evaluate_SVM(predicted, ref)
%--------------------------------------------------------------------------
% Modified by : Celine Provins (05.2020)
%
% Compute the matching ratio between the predicted label and the reference.
%
% INPUT
% predicted : vector containing the predicted labels
% ref : vector containing the true labels
%
% OUTPUT
% accuracy : float corresponding to the accuracy at predicting the right
%   labels
%--------------------------------------------------------------------------
    accuracy = 0.0;
    for i = 1:length(predicted)
        if strcmp(predicted(i),ref(i))
            accuracy = accuracy + 1.0;
        end
    end
    accuracy = accuracy /length(predicted);
end