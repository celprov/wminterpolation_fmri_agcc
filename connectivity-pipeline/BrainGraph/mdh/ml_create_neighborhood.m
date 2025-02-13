function [dirs,ndirs,norms,npar] = ml_create_neighborhood(dim,rempar)
% ML_CREATE_NEIGHBORHOOD Creates direction vectors for a 3D neighborhood.
%   [dirs,ndirs,norms] = ML_CREATE_NEIGHBORHOOD(dim)
%       dim specifies the size of the neighborhood cube, eg. dim=3 would
%       define a 26 neighborhood.
%   [dirs,ndirs,norms] = ML_CREATE_NEIGHBORHOOD(dim,rempar)
%       rempar specifies whether or not to remove parallel vectors, eg.
%       ML_CREATE_NEIGHBORHOOD(5,true) would return 98 and not 124 vectors.
%
%   Author:
%       Martin Larsson
%       June 2017

    r = (dim-1)/2;
    [y,x,z] = meshgrid(-r:r,-r:r,-r:r);
    not_origin = true(1,numel(x));
    not_origin((numel(x)+1)/2) = false;
    dirs = [x(not_origin); y(not_origin); z(not_origin)];

    if nargin > 1 && rempar || nargout > 1
        % Normalized vectors.
        norms = sqrt(sum(dirs.^2));
        ndirs = dirs./repmat(norms,size(dirs,1),1);
    end
    
    if nargout > 3 || nargin > 1 && rempar
        % Find parallel vectors.
        same = ndirs'*ndirs > 1-1e-6;
        include = false(size(dirs,2),1);
        for i=1:size(dirs,2)
            include(i) = norms(i) <= min(norms(same(i,:)));
        end
        npar = sum(include);
    end

    if nargin > 1 && rempar
        % Remove parallel vectors.
        dirs = dirs(:,include);
        ndirs = ndirs(:,include);
        norms = norms(include);
    end
end

