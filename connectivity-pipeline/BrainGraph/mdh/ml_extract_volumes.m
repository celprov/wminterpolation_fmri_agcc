function ml_extract_volumes(seg,varargin)
% ML_EXTRACT_VOLUMES Extracts NIFTI-volumes from a segmentation.
%   ML_EXTRACT_VOLUMES(seg,...)
%       seg is the file name of a segmentation NIFTI-file. Additional
%       arguments specify pairs of file names and ID-vectors. A NIFTI-file
%       will be saved to the specified location containing the voxels in
%       seg matching the IDs provided in the ID-vector.
%
%   Examples:
%       ML_EXTRACT_VOLUMES('aparc+aseg.nii','ribbon.nii',1000:2999);
%       ML_EXTRACT_VOLUMES('aparc+aseg.nii','thalamus.nii',[9 10 48 49]);
%
%   Author:
%       Martin Larsson
%       June 2017

    if nargin <= 1
        return;
    end

    [segh,segv] = ml_load_nifti(seg);

    h = segh;
    h.dt(1) = 2; % Data type: unsigned char (the smallest) 12/2017 D
    for i=1:2:length(varargin)
        h.fname = varargin{i};
        v = ismember(segv,varargin{i+1});
        spm_write_vol(h,v);
    end
end

