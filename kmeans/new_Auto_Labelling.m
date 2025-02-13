function atlasedVol = new_Auto_Labelling(segmentFiles, atlasFile, deform, outputDir, thresh)

% Syntax :
% Atlas =   new_Auto_Labelling(segmentFiles, atlasFile, deform, outputDir, thresh)
%
% This function returns an individual Atlas in Analyze format.
% The first step is to find the voxels that belong to gray matter tissue
% using the tissues probabilities maps previously obatained by the segmentation
% process. Due to the thresholding process, some holes, as well as isolated
% points are present at gray matter volume. To solve this problem is used the
% matlab function "imfill" to refills the internal holes and an internal
% function to reduce the isolated points(Uncoment if you want to employ it).
% Then each gray matter voxel is labeled with one structure label using the
% space transformation matrix (obtained in the normalization step), and an
% anatomical atlas (constructed by manual segmentation for a group of
% subjects).
%
% Input Parameters:
%  segmentFiles   : Cell array with names of the files containing the
%                   different brain segments. The first element of the
%                   array has to be the gray matter, for the other segments
%                   there is no prefered order, also the number of segments
%                   is not fix.
%  atlasFile      : Reference Atlas File.
%  deform         : Can be either a .mat deformation file, such as produced
%                   for example by SPM5/8 segmentation, or a .nii
%                   deformation field which can be produced by the SPM8 new
%                   segmentation.
%                   Note: Currently the corresponding algorithms for the 2
%                   files are implemented differently, if you use a mat
%                   file it has to be the mapping from the native space to
%                   the atlas. If you use a .nii deformation field, the
%                   mapping has to be from the atlas to the native space.
%  outputDir      : Output directory for segmented, normalized and atlased
%                   files. If the user doesn't change the output directory,
%                   the resulting files are saved in the same address than the
%                   Gray Matter Segmentation File.
%  thresh         : Threshold for gray matter segmentation.
%                   Just the voxels with  higher probability than the threshold
%                   are taken into acount in the automatic labelling step. If
%                   the threshold isn't specified then an automatic one is taken.
%                   All voxels with higher gray matter probabillity than 1-(GM+WM+CSF+...)
%                   are taken into account in the automatic labelling step.
%                   Being:
%                   GM(A voxel V belongs to gray matter tissue with a probability GM).
%                   WM(A voxel V belongs to white matter tissue with a probability WM).
%                   CSF(A voxel V belongs to cerebral spinal fluid with a probability CSF).
% Output Parameter:
%   atlasedVol   : Individual Gray Matter Atlas.
%
% Related References:
% 1.- Ashburner J, Friston K. Multimodal image coregistration and partitioning--
%     a unified framework. Neuroimage. 1997 Oct;6(3):209-17.
% 3.- Evans AC, Collins DL, Milner B (1992). An MRI-based Stereotactic Brain
%     Atlas from 300 Young Normal Subjects, in: Proceedings of the 22nd Symposium
%     of the Society for Neuroscience, Anaheim, 408.
%
%
% See also: imfill  imerode  imdilate  spm_segment atlasing spm_normalise
%__________________________________________________________________________
% Authors:  Yasser Aleman Gomez & Lester Melie Garcia
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2005
% Version $1.0
%__________________________________________________________________________
% 18.10.2010 Modified for spm8 new segmentation deformation fields(Manuel
%               Wuthrich)
% 20.10.2011 Fixed interpolation parameter which would produce weird regions
%           (Jonas Richiardi)
% 28.12.2011 forward and backwards compatibility of fileparts call (JR)
%=====================Checking Input Parameters===========================%

[pth,nam,ext] = spm_fileparts(segmentFiles{1});
%mkdir(outputDir);

mask = compMask(segmentFiles, thresh);

%-------- check if deform is a deformation field or a mat file ------
if length(deform) <= 4
    error(['The deformation file "' deform '" can not be read.']);
elseif strcmp('.mat', deform(end-3:end))
    atlasedVol = classicApply(atlasFile, deform, outputDir, mask);
elseif strcmp(ext, deform(end-3:end))
    atlasedVol = newApply(segmentFiles, atlasFile, deform, outputDir, mask);
else
    error(['The file extension "' deform(end-3:end) '" of the deformation file is unknown.']);
end
return;


function mask = compMask(segmentFiles, thresh)
%
% Input Parameters:
%   segmentFiles   : Cell array of filenames with different brain sagments,
%                    first has to be gray matter 
%   thresh         : Threshold, any voxel from the gray matter segment with
%                    a probability below thresh will be considered as not gray matter.
%
%  Output Parameters:
%   mask           : Mask with size of the segment files, 1 for gray matter
%                    and 0 for not gray matter.
%
%
% Note: This is based on the function 'Auto_Labelling' developed by 
%       Yasser Alem?n G?mez & Pedro Vald?s Hern?ndez.
%
%__________________________________________________________________________
% Author:  Manuel W?thrich
% EPFL
% Last update: December 20th 2010
% Version $1.0


%% ------ create array with all the segments ------------------------
segment = cell(size(segmentFiles,1),1);
total = 0; % sum of all probabilities in all TPMs
for i = 1:size(segmentFiles,1)
    segment{i} = spm_read_vols(spm_vol(segmentFiles{i}));
    if(i == 1), total = segment{1};
    else total = total + segment{i};
    end
end

%% -------------- create the mask for the gray matter ---------------
mask = true(size(segment{1}));
for i = 2:size(segmentFiles,1)
    mask(segment{1} <= segment{i}) = 0; % remove voxels where GM is less probable than other tissue classes
end

if  ~exist('thresh','var')||(isempty(thresh))||(thresh == 0)
    mask(segment{1} <= 1-total) = 0;
else
    mask(segment{1} < thresh) = 0;
end

%% ------- do morphology to fill holes and reduce isolated points --------
disp('Refilling...');
mask = imfill(logical(mask),'holes');
clear segment total segmentVol None;
mask = Iso_Rem(mask,7);
return

function atlasedVol = newApply(segmentFiles, atlasFile, deform, outputDir, mask)
%
% Input Parameters:
%  segmentFiles   : Cell array of filenames with different brain sagments,
%                   first has to be gray matter 
%  atlasFile      : Reference Atlas File.
%  deform         : Can be either a .mat deformation file, such as produced
%                   for example by SPM5/8 segmentation, or a .nii
%                   deformation field which can be produced by the SPM8 new
%                   segmentation.
%                   Note: Currently the corresponding algorithms for the 2
%                   files are implemented differently, if you use a mat
%                   file it has to be the mapping from the native space to
%                   the atlas. If you use a .nii deformation field, the
%                   mapping has to be from the atlas to the native space.
%  outputDir      : Output directory for segmented, normalized and atlased
%                   files. If the user doesn't change the output directory,
%                   the resulting files are saved in the same address than the
%                   Gray Matter Segmentation File.
%   mask           : Mask with size of the segment files, 1 for gray matter
%                    and 0 for not gray matter.
%
%  Output Parameters:
%   atlasedVol    : Individual Gray Matter Atlas.
%
%
%__________________________________________________________________________
% Author:  Manuel W?thrich
% EPFL
% Last update: Oct 20, 2011
% Version $1.1
% - fixed a critical bug with interpolation...

%% -------------- load deformation field -------------------------
[pth,~,ext] = spm_fileparts(segmentFiles{1});

defFieldFiles   = [repmat(deform,3,1), [',1,1';',1,2';',1,3']];
defFieldVol     = spm_vol(defFieldFiles);
deform          = cell(3,1);
deform{1}       = spm_load_float(defFieldVol(1));
deform{2}       = spm_load_float(defFieldVol(2));
deform{3}       = spm_load_float(defFieldVol(3));
mat             = defFieldVol(1).mat;

%intrp = [[1 1 1], 0 0 0];
% use 0th degree b-spline = NN interpolation in the SPM implementation.
% degree 1 == "derivatives of trilinear interpolation" which causes very
% weird regions to appear...
intrp = [[0 0 0], 0 0 0];

atlasVol = spm_vol(atlasFile);
M = inv(atlasVol.mat);

%% apply deformation field to atlas, mask the atlas and write it to file
% max-compatibility version of fileparts (old releases don't have ~)
if verLessThan('matlab', '7.13.0')
    [foo_atlasPath, foo_atlasBasename, foo_atlasExt, foo_atlasVer] = fileparts(atlasFile);
else
    [foo_atlasPath, foo_atlasBasename, foo_atlasExt] = fileparts(atlasFile); 
end
% nam = 'wNorm';
% ofnames = fullfile(outputDir,[nam,'_Atlas',foo_atlasBasename,ext]);
nam = 'wNorm_';
ofnames = fullfile(outputDir,[nam,foo_atlasBasename,ext]);
atlasedVol = struct('fname',ofnames,...
    'dim',[size(deform{1},1) size(deform{1},2) size(deform{1},3)],...
    'dt',atlasVol.dt,...
    'pinfo',atlasVol.pinfo,...
    'mat',mat,...
    'n',atlasVol.n,...
    'descrip',atlasVol.descrip);
C  = spm_bsplinc(atlasVol,intrp);
atlasedVol = spm_create_vol(atlasedVol);
d = cell(3,1);
for j=1:size(deform{1},3)
    d0    = {double(deform{1}(:,:,j)), double(deform{2}(:,:,j)),double(deform{3}(:,:,j))};
    d{1}  = M(1,1)*d0{1}+M(1,2)*d0{2}+M(1,3)*d0{3}+M(1,4);
    d{2}  = M(2,1)*d0{1}+M(2,2)*d0{2}+M(2,3)*d0{3}+M(2,4);
    d{3}  = M(3,1)*d0{1}+M(3,2)*d0{2}+M(3,3)*d0{3}+M(3,4);
    dat   = spm_bsplins(C,d{:},intrp);
    atlasedVol = spm_write_plane(atlasedVol,mask(:,:,j) .* dat,j);
end;
return;

%=====================Internal functions==================================%
function atlasedVol = classicApply(atlasFile, deform, outputDir, mask)
%
% Input Parameters:
%   atlasFile   : Reference Atlas File
%   matname     : Normalisation Transform File (from segmentation)
%   mask   : Mask for gray matter
%   outputDir  : Output Directory
%
% Note: This is based on the function 'spm_write_sn' developed by PhD.John Ashburner(FIL,UCL).
%
%__________________________________________________________________________
% Authors:  Yasser Alem?n G?mez & Pedro Vald?s Hern?ndez
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2005
% Version $1.0
warning off
defaults.analyze.flip = 0; % Not flipped Images
global defaults
atlasVol = spm_vol(atlasFile);
atlas = uint16(spm_read_vols(atlasVol));
if length(atlasVol.dim)==4
    dt = [atlasVol.dim(4) 0];
elseif length(atlasVol.dim)==3
    dt = atlasVol.dt;
end
% VAtlas_mat = atlasVol.mat;
load('-mat',deform);
[Vpth,Vname,Vext] = fileparts(VF.fname);
% JORI FIX FOR PC FILES ON MACS
bsidx=strfind(Vname,'\');
sidx=strfind(Vname,'/');
if ~isempty(bsidx) % we picked up PC file name paths
    Vname=Vname(bsidx(end)+1:end);
elseif ~isempty(sidx)
    error('code not ready to deal with macpaths on a PC');
end

mat = VF.mat;
if strcmp(spm('ver'),'SPM2')
    dim = [VF.dim(1:3) dt(1)];
elseif (strcmp(spm('ver'),'SPM5') || strcmp(spm('ver'),'SPM8'))     %adapted for spm8 by Manuel Wuthrich
    dim = [VF.dim(1:3)];
end
atlasedVol = struct('fname','','dim',dim,'mat',mat,'pinfo',[1 0 0]',...
    'descrip','Atlas image','dt',dt);
[foo_atlasPath, foo_atlasBasename, foo_atlasExt] = fileparts(atlasFile); 
atlasedVol.fname =[char(outputDir) filesep Vname '_Atlas' foo_atlasBasename Vext];
atlasedVol = spm_create_vol(atlasedVol);
x = 1:VG(1).dim(1); y = 1:VG(1).dim(2); z = 1:VG(1).dim(3); %dimensions of graymatter template
if ~isempty(Tr)
    BX = spm_dctmtx(VG(1).dim(1),size(Tr,1),x-1);
    BY = spm_dctmtx(VG(1).dim(2),size(Tr,2),y-1);
    BZ = spm_dctmtx(VG(1).dim(3),size(Tr,3),z-1);
end
[X,Y] = ndgrid(x,y); clear x y
y1 = single(0); y1(VG(1).dim(1),VG(1).dim(2),VG(1).dim(3)) = 0;
y2 = single(0); y2(VG(1).dim(1),VG(1).dim(2),VG(1).dim(3)) = 0;
y3 = single(0); y3(VG(1).dim(1),VG(1).dim(2),VG(1).dim(3)) = 0;
M = VG(1).mat;
for j=1:length(z);
    if ~isempty(Tr)
        X1 = X    + BX*get_2Dtrans(Tr(:,:,:,1),BZ,j)*BY';
        Y1 = Y    + BX*get_2Dtrans(Tr(:,:,:,2),BZ,j)*BY';
        Z1 = z(j) + BX*get_2Dtrans(Tr(:,:,:,3),BZ,j)*BY';
    else
        X1 = X; Y1 = Y; Z1 = z(j);
    end
    y1(:,:,j) = single(M(1,1)*X1 + M(1,2)*Y1 + M(1,3)*Z1 + M(1,4));
    y2(:,:,j) = single(M(2,1)*X1 + M(2,2)*Y1 + M(2,3)*Z1 + M(2,4));
    y3(:,:,j) = single(M(3,1)*X1 + M(3,2)*Y1 + M(3,3)*Z1 + M(3,4));
end
clear X1 Y1 Z1 X Y z
disp(['Inverting the deformations field...']);
M = Affine/VG(1).mat; M(4,:) = [0 0 0 1];
[iy1,iy2,iy3] = spm_invdef(y1,y2,y3,VF.dim(1:3),M,VG(1).mat);
clear y1 y2 y3
M = inv(atlasVol.mat);
for j = 1:VF.dim(3)
    A = zeros(VF.dim(1),VF.dim(2));
    disp(['Slice ----> ' num2str(j)]);
    if sum(sum(mask(:,:,j  )))~=0
        X2 = M(1,1)*double(iy1(:,:,j)) + M(1,2)*double(iy2(:,:,j)) + M(1,3)*double(iy3(:,:,j)) + M(1,4);
        Y2 = M(2,1)*double(iy1(:,:,j)) + M(2,2)*double(iy2(:,:,j)) + M(2,3)*double(iy3(:,:,j)) + M(2,4);
        Z2 = M(3,1)*double(iy1(:,:,j)) + M(3,2)*double(iy2(:,:,j)) + M(3,3)*double(iy3(:,:,j)) + M(3,4);
        A = spm_sample_vol(atlas,X2,Y2,Z2,0);
    end
    spm_write_plane(atlasedVol,squeeze(A.*mask(:,:,j)),j);
end
fclose all;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T2 = get_2Dtrans(T3,B,j)
d   = [size(T3) 1 1 1];
tmp = reshape(T3,d(1)*d(2),d(3));
T2  = reshape(tmp*B(j,:)',d(1),d(2));
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function I = Iso_Rem(T,Nhood)
%
%This function removes isolated points from white matter mask.
%
% Input Parameters:
%   T            : White Matter  Mask
%   Nhood        : Minimun number of neighbors.
% Output Parameters:
%   I            : White Matter Mask without isolated points
%__________________________________________________________________________
% Authors:  Yasser Alem?n G?mez
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2005
% Version $1.0

warning off
I = zeros(size(T)+2);
I(2:end-1,2:end-1,2:end-1) = T;
clear T
ind = find(I>0);
[x,y,z] = ind2sub(size(I), ind);
s = size(x,1);
sROI = zeros(size(I));
% figure;spy(I(:,:,160));
[X, Y, Z] = meshgrid(-1:1,-1:1,-1:1);
X = X(:);Y = Y(:);Z = Z(:);
Neib = [X Y Z];clear X Y Z;
pos = find((Neib(:,1)==0)&(Neib(:,2)==0)&(Neib(:,3)==0));
Neib(pos,:) = [];
for i =1:26
    M = Neib(i,:);
    S = [x y z] + M(ones(s,1),:);
    ind2 = sub2ind(size(I),S(:,1),S(:,2),S(:,3));
    sROI(ind) = sROI(ind) + I(ind2);
end
ind = find(sROI<Nhood);
I(ind) =0;
I = I(2:end-1,2:end-1,2:end-1);
return;


