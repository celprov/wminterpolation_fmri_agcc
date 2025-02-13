function [fileOutput] = da_corregister_reslice( ref, source, interp, prefix )
%DA_RESIZE Reslice a volume using SPM
%   ref: path of reference volume
%   source: path of volume to reslice
%   interp: 0 for nearest neighbor, 1 for trilinear interpolation
%   prefix: prefix added to the name of the resliced file

if nargin < 3
    interp = 0;
end
if nargin < 4
    prefix = 'r';
end

if interp ~=0 && interp ~=1
    error('Interpolation must be either 0 or 1')
end

job.ref = {strcat(ref,',1')};
job.source = {strcat(source,',1')};
job.roptions.interp = interp;
job.roptions.wrap = [0 0 0];
job.roptions.mask = 0;
job.roptions.prefix = prefix;

out = spm_run_coreg(job);
fileOutput = strsplit(out.rfiles{1},',');
fileOutput = fileOutput{1};
end

