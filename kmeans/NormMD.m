function [out,idx_nonzero] = NormMD(isBC,idx_active, V, defld,...
    indices_wb,hdr,fHeader,subInpaintingPath)
%--------------------------------------------------------------------------
% Modified by : Celine Provins (04.2020)
%
% INPUT
% isBC : bool determining whether the data are BC or GSR
% idx_active : vector containing the frame indices that we consider active
% V : matrix containing the volumes (BC or GSR)
% defld : path to the deformation field (string)
% indices_wb : non-zero indices in the inpainted volumes
% hdr : string corresponding to the header of the model nifti which
%   determines the dimension of the saved nifit
% fHeader : matrix determining the dimension of the saved nifti
% subInpaintingPath : string corresponding to the path towards the folder
%   containing the inpainted folder for one particular subject
%
% OUTPUT
% out : matrix containing the normalized active frames 
% idx_nonzero : vector containing the linear coordinates of the non-zero
%   voxels in the normalized frames
%--------------------------------------------------------------------------
    for i = 1:length(idx_active)
        %Write the selected frames one by one
        A=zeros(fHeader.dim);
        A(indices_wb) = V(idx_active(i),:);
        if isBC
            if ~exist(fullfile(subInpaintingPath,strcat('BCsel',num2str(i,'%03d'),'.nii')),'file')         
                cbiWriteNifti(fullfile(subInpaintingPath,...
                    strcat('BCsel',num2str(i,'%03d'),'.nii')),A,hdr,'float32');
            end
        else
            if ~exist(fullfile(subInpaintingPath,strcat('INsel',num2str(i,'%03d'),'.nii')),'file')       
                cbiWriteNifti(fullfile(subInpaintingPath,...
                strcat('INsel',num2str(i,'%03d'),'.nii')),A,hdr,'float32');
            end
        end
    end

    %Normalize the selected frames to MNI space
    if isBC
        if ~exist(fullfile(subInpaintingPath,strcat('wBCsel',num2str(i,'%03d'),'.nii')),'file')
            NormalizeToMD(subInpaintingPath,'BCsel',defld);
        end
    else
        if ~exist(fullfile(subInpaintingPath,strcat('wINsel',num2str(i,'%03d'),'.nii')),'file')
            NormalizeToMD(subInpaintingPath,'INsel',defld);
        end
    end

    %Read the normalized data out
    for i = 1:length(idx_active)
        if isBC
            A = spm_read_vols(spm_vol(fullfile(subInpaintingPath,...
                strcat('wBCsel',num2str(i,'%03d'),'.nii'))));
            out(i,:)= A(:);
        else
            A = spm_read_vols(spm_vol(fullfile(subInpaintingPath,...
                strcat('wINsel',num2str(i,'%03d'),'.nii'))));
            out(i,:)= A(:);
        end
    end
    out(isnan(out)==1)=0;
    idx_nonzero = find(out(1,:)~=0);
end
