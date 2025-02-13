function [mat_file_name] = ml_fibgz2mat(file_name)
% ML_FIBGZ2MAT Converts FIB.GZ-files into MAT-files.
%   ML_FIBGZ2MAT(file_name)
%       file_name specifies the FIB.GZ-file to convert. If no file name is
%       provided a file dialog will be opened. The name of the new file is
%       returned.
%
%   Author:
%       Martin Larsson
%       May 2017

    if ~exist('file_name', 'var') || isempty(file_name)
        [file_name, path_name] = uigetfile({'*.fib.gz', 'Compressed FIB-files'});
        if file_name == 0
            return;
        end
        file_name = fullfile(path_name, file_name);
    end
    
    [pathstr,name,~] = fileparts(file_name);
    fib_file_name = fullfile(pathstr, name);
    [pathstr,name,~] = fileparts(fib_file_name);
    mat_file_name = fullfile(pathstr, strcat(name, '.mat'));
    gunzip(file_name);
    movefile(fib_file_name, mat_file_name);
end

