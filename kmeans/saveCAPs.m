function [] = saveCAPs(prefix,CAP,hdr,k,idx_nonzero,isBC,idx_col,capPath)
%--------------------------------------------------------------------------
% Modified by : Céline Provins (04.2020)
% 
% Save the CAPs in nifti format
%
% prefix : string precising whether its BC or GSR CAP
% CAP : matrix containing the CAPs to save
% hdr : nifti model's header determining the dimension of the CAPs to save
% k : number of clusters
% idx_nonzero : vector containing the non-zero voxel linear coordinate
% isBC : bool precising whether we work with BC data (1) or GSR (0)
% idx_col : vector containing the indices matching between CAP_BC and
%   CAP_GSR
% capPath : string containing the path where to save the CAPs
%--------------------------------------------------------------------------

V = zeros(hdr.dim(2),hdr.dim(3),hdr.dim(4));
for l = 1:k
    if isBC
        V(idx_nonzero) = zscore(CAP(l,:));
    else 
        V(idx_nonzero) = zscore(CAP(idx_col(l),:));
    end
    cbiWriteNifti(fullfile(capPath,strcat(prefix,'CAP_',num2str(l),'_mean.nii')),V,hdr,'float32');
end
