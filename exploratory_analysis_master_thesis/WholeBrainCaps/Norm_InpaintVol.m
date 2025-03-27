function [out,out_lowRes,idx_nonzero,idx_nonzero_lowRes] = ...
Norm_InpaintVol(isBC, prefix, deformationPrefix,BrainGraphPath,...
InpaintingPath,subject,segPath)
%--------------------------------------------------------------------------
% Created by : Celine Provins (06.2020)
%
% Write all the volumes one by one separately and normalize those volumes
% in MNI space
%
% INPUT
% isBC : bool that precise whether you work with BC or GSR data
% prefix : prefix of the file where the volumes are stored
% deformatinPrefix : string containing the prefix of the deformation field
% BrainGraphPath : path towards the folder containing the brain graphs
% InpaintingPath : path towards the folder containing the inpainted volumes
% subject : name of the subject 
% segPath : path towards the folder containing the segmented files of the
%   subject
%
% OUTPUT
% out : matrix containing the normalized active frames
% out_lowRes : matrix containing the normalized active frames in low
%   resolution
% idx_nonzero : vector containing the linear coordinates of the non-zero
%   voxels in the normalized frames
% idx_nonzero_lowRes : vector containing the linear coordinates of the non-zero
%   voxels in the normalized frames in low resolution
%--------------------------------------------------------------------------
    %% Load variables
    indices_wb = load(fullfile(BrainGraphPath,subject,'ODF_Neigh_3_ODFPower_40','Spectrum_WB_improved','indices_wb.mat'));
    indices_wb = indices_wb.indices_wb; %non-zero indices of the inpainted brain
    
    subInpaintingPath = fullfile(InpaintingPath,subject);
    cd(subInpaintingPath);
    file = dir(strcat(prefix,'*.mat'));
    V = load(file.name);
    V = V.V;
    
    fHeader = load(fullfile(subInpaintingPath,'fHeader.mat'));
    fHeader = fHeader.fHeader;
    
    cd(segPath);
    defld = dir(deformationPrefix);
    deformfield = fullfile(defld.folder,defld.name);
  
    hdr = cbiReadNiftiHeader(fullfile(BrainGraphPath,subject,'T1w','new_brainmask.nii'));

    %% Detrend the volumes if hasn't been done before
    if isBC
        disp('Let us detrend volume..')
        CUTNUMBER=10;
        SegmentLength = ceil(size(V,2) / CUTNUMBER);
        for iCut=1:CUTNUMBER
            if iCut~=CUTNUMBER
                Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
            else
                Segment = (iCut-1)*SegmentLength+1 : size(V,2);
            end
            V(:,Segment) = detrend(V(:,Segment));
            V(:,Segment) = detrend(V(:,Segment),'constant');
            fprintf('.');
        end
        V= zscore(V);
    end

    %% Write the volumes
    for i = 1:size(V,1)
        A=zeros(fHeader.dim);
        A(indices_wb) = V(i,:);
        if isBC
            if ~exist(fullfile(subInpaintingPath,strcat('WholeBrainBC',num2str(i,'%03d'),'.nii')),'file')         
                cbiWriteNifti(fullfile(subInpaintingPath,...
                    strcat('WholeBrainBC',num2str(i,'%03d'),'.nii')),A,hdr,'float32');
            end
        else
            if ~exist(fullfile(subInpaintingPath,strcat('WholeBrainIN',num2str(i,'%03d'),'.nii')),'file')       
                cbiWriteNifti(fullfile(subInpaintingPath,...
                strcat('WholeBrainIN',num2str(i,'%03d'),'.nii')),A,hdr,'float32');
            end
        end
    end
    
    %% Normalize the volumes
    if isBC
        if ~exist(fullfile(subInpaintingPath,strcat('w3WholeBrainBC',num2str(i,'%03d'),'.nii')),'file')
            NormalizeToMNI_size(subInpaintingPath,'WholeBrainBC',deformfield,3);
        end
        if ~exist(fullfile(subInpaintingPath,strcat('w2WholeBrainBC',num2str(i,'%03d'),'.nii')),'file')
            NormalizeToMNI_size(subInpaintingPath,'WholeBrainBC',deformfield,2);
        end
    else
        if ~exist(fullfile(subInpaintingPath,strcat('w3WholeBrainIN',num2str(i,'%03d'),'.nii')),'file')
            NormalizeToMNI_size(subInpaintingPath,'WholeBrainIN',deformfield,3);
        end
        if ~exist(fullfile(subInpaintingPath,strcat('w2WholeBrainIN',num2str(i,'%03d'),'.nii')),'file')
            NormalizeToMNI_size(subInpaintingPath,'WholeBrainIN',deformfield,2);
        end
    end

    %% Read the normalized data out
    for i = 1:size(V,1)
        if isBC
            A = spm_read_vols(spm_vol(fullfile(subInpaintingPath,...
                strcat('w3WholeBrainBC',num2str(i,'%03d'),'.nii'))));
            B = spm_read_vols(spm_vol(fullfile(subInpaintingPath,...
                strcat('w2WholeBrainBC',num2str(i,'%03d'),'.nii'))));
            out_lowRes(i,:)= A(:);
            out(i,:)= B(:);
        else
            A = spm_read_vols(spm_vol(fullfile(subInpaintingPath,...
                strcat('w3WholeBrainIN',num2str(i,'%03d'),'.nii'))));
            B = spm_read_vols(spm_vol(fullfile(subInpaintingPath,...
                strcat('w2WholeBrainIN',num2str(i,'%03d'),'.nii'))));
            out_lowRes(i,:)= A(:);
            out(i,:)= B(:);
        end
    end
    out(isnan(out)==1)=0;
    out_lowRes(isnan(out_lowRes)==1)=0;
    idx_nonzero = find(out(1,:)~=0);
    idx_nonzero_lowRes = find(out_lowRes(1,:)~=0);
end
