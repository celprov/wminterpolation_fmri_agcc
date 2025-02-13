function [ fileOutput ] = da_resize( fileInput, fileRef )
%DA_RESIZE_AND_CLEAN Resize a nifti file to match the size of a reference
%file
%   This function resizes an input nifti file to match the size of a
%   reference file. 
%
%   Inputs:
%       fileInput: full address of nifti file to resize.
%       fileRef: full address of nifti file used as reference for size and
%           coregistration.
%   
%   Outputs:
%       outputFile: full address of created file, with naming convention
%       inputFileName_referenceVoxelSize.nii, and created in the same
%       folder as the input file.
%
%   David Abramian, May 2017

% Check input for errors
if ~exist(fileInput,'file')
    error('Input file does not exist')
end

if ~exist(fileRef,'file')
    error('Reference file does not exist')
end

[dirInput, name, ext] = fileparts(fileInput);
if ~strcmp(ext, '.nii')
    error('Input must be a nifti file')
end

[~, ~, ext] = fileparts(fileRef);
if ~strcmp(ext, '.nii')
    error('Reference must be a nifti file')
end


% Determine output voxel size
hRef = spm_vol(fileRef);
outputVoxelSize = abs(hRef(1).mat(1,1));

if rem(outputVoxelSize,1) == 0
    strOutputVoxelSize = num2str(outputVoxelSize,'%.1f');
else
    strOutputVoxelSize = num2str(outputVoxelSize);
end

% Name of output file
fileOutput = fullfile(dirInput, [name,'_', strOutputVoxelSize, ext]);

% Resize file
da_corregister_reslice(fileRef, fileInput, 0, 'rnn_');
movefile(fullfile(dirInput, ['rnn_', name, ext]), fileOutput); 


