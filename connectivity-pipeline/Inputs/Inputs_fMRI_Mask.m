
% Paths
param.SaveInpainting = fullfile(param.HCPDatapath2, 'Inpainting_results_test2');

Mask = spm_read_vols(spm_vol(fullfile(param.structseg,...
    ['c1', param.T1file])));

fHeader = spm_vol(fullfile(param.structseg,...
    ['c1', param.T1file]));
param.fHeader = fHeader;

hdr=cbiReadNiftiHeader(fHeader.fname);

Mask(Mask<0.3)=0;

Mask = bwareaopen(logical(Mask),30);

Mask = imfill(Mask,'holes');

BCBrainMask = logical(Mask);

param.smoothR = 6;

indices_wb = param.ODF.G.indices_wb;

BC = find(BCBrainMask);

[sim_ind,~,BC] = intersect(BC, indices_wb);
