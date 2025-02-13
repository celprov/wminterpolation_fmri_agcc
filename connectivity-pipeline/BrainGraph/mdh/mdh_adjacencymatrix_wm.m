function [A,indices,G,O] = mdh_adjacencymatrix_wm(G,odfPow,neighb,dsi_root,varargin)
% MDH_ADJACENCYMATRIX_WM Constructs white matter graph A matrix.
%
%
% 
% 
% Contributers: 
%    Martin Larsson
%    David Abramian
%    Hamid Behjat
% 2017-2018. 

cntr = {
    'N_fibers',5,...     % maximum count of the resolving fibers
    };
argselectAssign(cntr);
argselectCheck(cntr,varargin);
argselectAssign(varargin);

fprintf('\n.Constructing WM graph adjacency matrix..\n')

f = G.f;

[wHeader,wVol] = ml_load_nifti(f.white_125);

% Build ODFs. 
f = mdh_dsistudio(f,dsi_root,'N_fibers',N_fibers);

% Build A. 
[A,odf,index,odf_vertices,odf_faces,f] = mdh_adjacencymatrix_wm_pvt(...
    f.odfs,wVol,neighb,odfPow,f);

indices = find(wVol);
G.N_wm = size(A,1);
G.A_wm = A;
G.A_wm_uniform = double(logical(A)); % non-weighted A
G.indices_wm = indices;
G.dim = wHeader.dim;
[G.xyz(:,1), G.xyz(:,2), G.xyz(:,3)] = ind2sub(G.dim, indices);
G.f = f;

O = struct;
O.odf_vertices = odf_vertices';
O.odf_faces = odf_faces';
O.index = index(:,indices);
O.odf = odf;

A_wm = A; %#ok<NASGU>
indices_wm = indices; %#ok<NASGU>
save(fullfile(G.f.graph,'A_wm.mat'),'A_wm');
save(fullfile(G.f.graph,'indices_wm.mat'),'indices_wm');
fprintf('.done.\n')
end

function [A,odf,index,odf_vertices,odf_faces,f] = ...
    mdh_adjacencymatrix_wm_pvt(f_odf,mask,ndim,odf_power,f)

dim = size(mask);

% Convert FIB to MAT.
assert(strcmp(f_odf(end-6:end),'.fib.gz'),'Expected FIB.GZ-file for ODFs.');
f_odf = f_odf(1:end-length('.fib.gz'));
f_odf_fib = [f_odf '.fib.gz'];
f_odf_mat = [f_odf '.mat'];

if ~exist(f_odf_mat,'file')
    ml_fibgz2mat(f_odf_fib);
end

% Load ODFs.
[fa,odf,index,odf_vertices,odf_faces] = ml_dsi_load_odfs(f_odf_mat,true);
odf = odf.^odf_power;
odf = [odf; odf];

% Error check, misc.
if nnz(mask) ~= size(odf,2) || size(odf,1) ~= size(odf_vertices,2)
    error('Mismatching sizes.');
end
if ~isequal(mask(:)',fa(1,:)~=0)
    error('White matter masks misaligned.');
end

maski = find(mask);
voxels = length(maski);
maski_inv = zeros(dim);
maski_inv(maski) = 1:voxels;

% Find neighboring voxels and direction of edges.
[dirs,ndirs,~,npar] = ml_create_neighborhood(ndim);
neigh = size(dirs,2);

ci = repelem(maski,neigh,1);
cs = zeros(size(ci,1),3);
[cs(:,1),cs(:,2),cs(:,3)] = ind2sub(dim,ci);
ns = cs+repmat(dirs',voxels,1);
ni = sub2ind(dim,ns(:,1),ns(:,2),ns(:,3));
ni = reshape(ni,neigh,voxels);
ci = maski_inv(ci);
ni = maski_inv(ni);
lni = logical(ni);

% Construct graph.
%omega = 4*pi/neighborhood;
omega = 4*pi/npar; % We still want 98 for 5x5x5 neighborhood.
% 2*theta is a the apex angle of a cone.
costheta = 1-omega/(2*pi);
cone_set = ndirs'*odf_vertices > costheta;
cone_set = cone_set./repmat(sum(cone_set,2),1,size(cone_set,2));

% deltaS has no impact due to the normalization that follows.
%     deltaS = 4*pi/size(odf_vertices,2);
%     Pdiff = deltaS*cone_set*odf.*lni;

Pdiff = cone_set*odf.*lni;
Pdiff = 0.5*Pdiff./repmat(max(Pdiff),neigh,1);

A = sparse(ci(lni),ni(lni),Pdiff(lni),voxels,voxels);
A = A+A';

f.odf_fib = f_odf_fib;
f.odf_mat = f_odf_mat;
end
