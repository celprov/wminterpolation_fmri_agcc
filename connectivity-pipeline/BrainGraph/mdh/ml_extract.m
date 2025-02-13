function ml_extract(fname)
% ML_EXTRACT Extracts a file.
%   ML_EXTRACT(fname)
%       If the given files does not exists this function will look for
%       [fname '.gz'] and extract it.
%
%   Author:
%       Martin Larsson
%       June 2017

    if ~exist(fname,'file')
        comp_fname = [fname '.gz'];
        if ~exist(comp_fname,'file')
            error(['File does not exist: ' fname]);
        end
        gunzip(comp_fname);
    end
end

