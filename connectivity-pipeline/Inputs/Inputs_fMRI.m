%%%%%% INPUT FUNCTIONAL VOLUMES PATH


% filepath of preprocessed functional volumes:
% should already be in matrix form (timepoints x number of voxel)

param.functional = fullfile(param.HCPDatapath2, param.subject,'func','realigned');

prefix = 'Cov';

% to smooth volumes put 1
Smooth = 1;

param.smoothR = 6;

% to detrend volumes (put 1 if volumes were not detrended in preprocessing)
Detrend  = 1;

% Contains all inputs regarding functional MRI volumes which are subject-specific
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Functional volumes 

if ~exist(fullfile(param.SaveInpainting,['BCVolumes_',prefix,'_',param.subject,'.mat']),'file')

directory = dir(fullfile(param.functional, [prefix,'*']));


if Smooth 

    directory2 = dir(fullfile(param.functional, ['s',num2str(param.smoothR),prefix,'*']));
%         for i = 1:length(directory2)
%             delete(directory2(i).name)
%         end

    %directory2 = struct();
    l = length(directory2);
    if l == 0
        l = 1;
    end
    fVolsFNlist = struct2cell(directory); 
    fVolsFNlist = fVolsFNlist(1,:);
    fVolsFNlist_fullpath=cellfun(@(x) fullfile(param.functional,x), fVolsFNlist,'UniformOutput',false);
    fVolsFNlist_fullpath = fVolsFNlist_fullpath(l:end);
    disp(['There are ..',num2str(l),' already smoothed'])

    % Create batch
    spm_jobman('initcfg');
    matlabbatch{1}.spm.spatial.smooth.data = fVolsFNlist_fullpath';
    matlabbatch{1}.spm.spatial.smooth.fwhm = [param.smoothR param.smoothR param.smoothR];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = ['s',num2str(param.smoothR)];
    spm_jobman('run',matlabbatch);


end

directory = dir(fullfile(param.functional, ['s',num2str(param.smoothR),prefix,'*']));

V = zeros(length(directory), length(param.indices_wb));
for i = 1:length(directory)
    if ~mod(i,100) 
        disp(['Reading..',num2str(i),'th'])
    end
    temp = spm_read_vols(spm_vol(fullfile(param.functional,directory(i).name)));
    V(i,:) = temp(param.indices_wb);
end

save(fullfile(param.SaveInpainting,['BCVolumes_',prefix,'_',param.subject,'.mat']),'V','-v7.3')

else
    load(fullfile(param.SaveInpainting,['BCVolumes_',prefix,'_',param.subject,'.mat']))
end

% If the functional volumes were not detrended, we perform detrending
if Detrend
    disp('Let us detrend volume first..')
    CUTNUMBER=10;
    SegmentLength = ceil(size(V,2) / CUTNUMBER);
    for iCut=1:CUTNUMBER
        if iCut~=CUTNUMBER
            Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
        else
            Segment = (iCut-1)*SegmentLength+1 : size(V,2);
        end
        V(:,Segment) = detrend(V(:,Segment));
        V(:,Segment) = detrend(V(:,Segment),'constant');
        fprintf('.');
    end
end
% 
NumScans = size(V,1);
% Vz = zscore(V);

%% Print detrended volumes
% cbiPath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine/mrTools';
% addpath(genpath(cbiPath));
% 
% directory = dir(fullfile(param.functional, ['s',num2str(param.smoothR),prefix,'*']));
% for iter = 1:length(directory)  
%     temp = zeros(param.fHeader.dim);
%     temp(param.indices_wb) = V(iter,:);
%     nfile = fullfile(param.SaveInpainting,['DetrendedBefore_Volumes',sprintf('_%.3d',iter),'.nii']);
%     hdr=cbiReadNiftiHeader('/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine/data/RestingState/s012/func/realigned/s6Cov_regrfR649005-0019-00056-000056-01.nii');
%     cbiWriteNifti(nfile,temp,hdr,'float32');
% end