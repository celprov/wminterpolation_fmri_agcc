function [A,indices,G,O] = mdh_adjacencymatrix_wb2(G,odfPow,neighb,dsi_root,param,varargin)
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

[wHeader,wVol] = ml_load_nifti(f.brainmask_fs_125);

% Build ODFs. 
f = mdh_dsistudio_wb(f,dsi_root,'N_fibers',N_fibers);

% Build A. 
[A,odf,index,odf_vertices,odf_faces,f,mask] = mdh_adjacencymatrix_wb_pvt(...
    f.odfs,wVol,neighb,odfPow,f,param);

wVol = mask;
indices = find(wVol);
G.N_wb = size(A,1);
G.A_wb = A;
G.A_wb_uniform = double(logical(A)); % non-weighted A
G.indices_wb = indices;
G.dim = wHeader.dim;
[G.xyz(:,1), G.xyz(:,2), G.xyz(:,3)] = ind2sub(G.dim, indices);
G.f = f;

O = struct;
O.odf_vertices = odf_vertices';
O.odf_faces = odf_faces';
O.index = index(:,indices);
O.odf = odf;

A_wb = A; %#ok<NASGU>
indices_wb = indices; %#ok<NASGU>
save(fullfile(G.f.graph,'A_wb.mat'),'A_wb');
save(fullfile(G.f.graph,'indices_wb.mat'),'indices_wb');
fprintf('.done.\n')
end

function [A,odf,index,odf_vertices,odf_faces,f,mask] = ...
    mdh_adjacencymatrix_wb_pvt(f_odf,mask,ndim,odf_power,f,param)

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
[fa,odf,index,odf_vertices,odf_faces] = ml_dsi_load_odfs_wb(f_odf_mat,true);
odf = odf.^odf_power;
odf = [odf; odf];


% Quantitative anisotropy, consider only major fiber
Qa = max(fa,[],1);
Qa = Qa./max(Qa);
Qa(isnan(Qa)) = 0;

% Error check, misc.
if nnz(mask) ~= size(odf,2) || size(odf,1) ~= size(odf_vertices,2)
    disp('Mismatching sizes.');
    mask = reshape(Qa, dim);
    mask = logical(mask);
end
if ~isequal(mask(:)',fa(1,:)~=0)
    error('White matter masks misaligned.');
end


% check if there are voxels with no existing fiber (0 sum odf)
zerovox = sum(odf,1);
if(nnz(zerovox)~=size(odf,2))
    disp(['removing voxels with no odf values.. There are ..',...
        num2str(length(find(zerovox==0))),'.. of voxels to be removed'])
    odf(:,zerovox==0)=[];
    mask2 = zeros(dim);
    mask2(mask) = zerovox;
    mask = logical(mask2);
end


% Find neighboring voxels and direction of edges.
[dirs,ndirs,~,npar] = ml_create_neighborhood(ndim);
neigh = size(dirs,2);

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


maski = find(mask);
voxels = length(maski);
maski_inv = zeros(dim);
maski_inv(maski) = 1:voxels;

ci = repelem(maski,neigh,1);
cs = zeros(size(ci,1),3);
[cs(:,1),cs(:,2),cs(:,3)] = ind2sub(dim,ci);
ns = cs+repmat(dirs',voxels,1);
ni = sub2ind(dim,ns(:,1),ns(:,2),ns(:,3));
ni = reshape(ni,neigh,voxels);
ci = maski_inv(ci);
ni = maski_inv(ni);
lni = logical(ni);


%Save copy of QA and updated brainmask
h = spm_vol(f.white_125);
hdr = cbiReadNiftiHeader(h.fname);
niftifile = fullfile(f.extracted, 'QA_max.nii');
% if ~exist(niftifile,'file')   
    cbiWriteNifti(niftifile,reshape(Qa,dim),hdr,'float32')
% end
niftifile2 = fullfile(f.extracted, 'new_brainmask.nii');
% if ~exist(niftifile2,'file')
    cbiWriteNifti(niftifile2,mask,hdr,'float32')
% end

Pdiff = cone_set*odf.*lni;
Pdiff = 0.5*Pdiff./repmat(max(Pdiff),neigh,1);

Qa = Qa(mask);
Qa = repmat(Qa, [size(Pdiff,1), size(Qa,1)]);

% Note: Some voxels return nan values upon normalizing the weights locally
% (Anjali, 4 March 2018), putting them to 0.
% For some subjects, Pdiff gives 0 weights in some voxels, odf appears 0 (Anjali,5 March 2018)..
Pdiff(isnan(Pdiff)) = 0;

% Multiply weights with Qa, gamma filtered with alpha, takes only the major fiber
Pdiff = (Qa.^param.alpha).*Pdiff;

A = sparse(ci(lni),ni(lni),Pdiff(lni),voxels,voxels);
A = A+A';

% Check if graph is connected
conn = any(A);
if length(conn)==size(A,1)
    disp('Fully connected graph..')
else
    disp(['There are ..',num2str(length(conn)),'.. unconnected voxels'])
end
    

% Check for outlier voxels
% A2= A(find(A~=0));
% mean_allvoxel = mean(A2(:));
% small_num = eps; % order of magnitude
% thresh = mean_allvoxel*small_num;
% out = find(A2<thresh);


f.odf_fib = f_odf_fib;
f.odf_mat = f_odf_mat;
end
