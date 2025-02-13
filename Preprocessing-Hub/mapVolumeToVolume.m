function Vout=mapVolumeToVolume(Vin_fn,Vout_fn)
% map voxels in the space of input volume to voxels in the space
% of output volumes
% Vout_fn: uses only header info
%
% v1.0 Jonas Richiardi
% - initial release, based on code by Dimitri Van De Ville and Jonas
% Richiardi

% Read output space header
Vout_i=spm_vol(Vout_fn);
Vout=zeros(Vout_i.dim);

% read and load input space header and data
Vin_i=spm_vol(Vin_fn);
Vin=spm_read_vols(Vin_i);

% generate all coordinates in output space
[x1,x2,x3]=ndgrid(1:Vout_i.dim(1),1:Vout_i.dim(2),1:Vout_i.dim(3));
idx=1:numel(Vout); % map all voxels

% take every voxel in the volume spanned by the output images,
% compute its real-world position in mm, then map input image

oobList=zeros(0,4); % list of out-of-bound input voxels
for iter=1:length(idx),
    oob=false;
    % recover world-space position of this voxel in mm from affine
    % transform matrix
    mm=Vout_i.mat*[x1(idx(iter)) x2(idx(iter)) x3(idx(iter)) 1]';
    % convert this position into index of the closest structural voxel
    vx=round(Vin_i.mat\[mm(1) mm(2) mm(3) 1]'); 
    vx(vx<=0)=1;
    vxOri=vx;
    % remap out-of-bounds voxels to last  
    if vx(1)>Vin_i.dim(1), vx(1)=Vin_i.dim(1); oob=true; end 
    if vx(2)>Vin_i.dim(2), vx(2)=Vin_i.dim(2); oob=true; end
    if vx(3)>Vin_i.dim(3), vx(3)=Vin_i.dim(3); oob=true; end
    if (oob==true), oobList(end+1,:)=vxOri; end
    % idx(iter): current voxel
    Vout(idx(iter))=Vin(vx(1),vx(2),vx(3));
%     if any(Vout(idx(iter))<0)  %Vout
%         warning('mapV2V:negativeVal',['Negative voxel values at ' num2str(iter)]);
%     end
end;