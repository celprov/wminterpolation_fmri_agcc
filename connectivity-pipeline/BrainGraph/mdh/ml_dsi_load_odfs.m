function [fa,odf,index,odf_vertices,odf_faces] = ml_dsi_load_odfs(file_name,flipY)
% ML_DSI_LOAD_ODFS Load ODF data outputted from DSI Studio.
%   [fa,odf,index,odf_vertices,odf_faces] = ML_DSI_LOAD_ODFS(file_name)
%       file_name specifies the MAT-file that contains the ODF data. DSI
%       Studio outputs a FIB.GZ-file that can be converted to a MAT-file
%       using ml_fibgz2mat.
%   [fa,odf,index,odf_vertices,odf_faces] = ML_DSI_LOAD_ODFS(file_name,flipY,dim)
%       flipY specifies whether or not to flip the data along the Y-axis
%       (second dimention). This is required to match the voxel space of
%       DSI Studio (-X,-Y,Z) with the HCP dataset (-X,Y,Z).
%
%   See also ML_FIBGZ2MAT.
%
%   Author:
%       Martin Larsson
%       May 2017

    if ~exist('file_name', 'var') || isempty(file_name)
        [file_name, path_name] = uigetfile({'*.mat','MATLAB files'});
        if file_name == 0
            return;
        end
        file_name = fullfile(path_name, file_name);
    end

    load(file_name);

    ml_catvars('fa',1);
    ml_catvars('odf',2);
    ml_catvars('index',1);
    index = index+1;

    if nargin > 1 && flipY
        fai = find(fa(1,:)~=0);
        fai_inv = zeros(dimension);
        fai_inv(fai) = 1:length(fai);
        fai_inv = flip(fai_inv,2);
        perm = fai_inv(fai_inv~=0);
        odf = odf(:,perm);
        odf_vertices(2,:) = -odf_vertices(2,:);

        for i=1:size(fa,1)
            fa(i,:) = reshape(flip(reshape(fa(i,:),dimension),2),1,size(fa,2));
            index(i,:) = reshape(flip(reshape(index(i,:),dimension),2),1,size(index,2));
        end
    end
    
  
    
end

