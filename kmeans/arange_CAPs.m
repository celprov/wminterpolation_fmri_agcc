function [matrix_out2,idx_col] = arange_CAPs(matrix_in1,matrix_in2)
%------------------------------------------------------------------------
% Created by : Zhiwei Huang (10.2019)
%
% Match the indices of the CAP contained in matrix_in2 with the indices of
% the CAP in matrix_in1, so that similar CAPs carry the same number
%
% INPUT
% matrix_in1 : matrix containing the reference CAP
% matrix_in2 : CAP for which the indices are matched to matrix_in1
%
% OUTPUT
% matrix_out2 : matrix_in2 rearranged so that the CAP number match with
%   matrix_in1
% idx_col : vector containing the indices matching between CAP_BC and
%   CAP_GSR
%------------------------------------------------------------------------
similarity = zeros(size(matrix_in1,1),size(matrix_in2,1));
for i = 1:size(matrix_in1,1)
    for j = 1:size(matrix_in2,1)
        similarity(i,j) = pdist2(matrix_in1(i,:),matrix_in2(j,:),'cosine');
    end
end

% find the assignment
[idx_col,~] = munkres(similarity);
matrix_out1 = matrix_in1;
matrix_out2 = zeros(size(matrix_in2));
for i = 1:size(matrix_in1,1)
    matrix_out2(i,:) = matrix_in2(idx_col(i),:);
end