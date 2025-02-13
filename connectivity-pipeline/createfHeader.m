function createfHeader(dataset)
    HCPDatapath = '/media/miplab-nas2/Data/AgCC/Anjali/BrainGraph_results';
    InpaintingPath = fullfile('/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine/data', dataset,'Inpainting_results');
    spmPath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine/spm12';
    addpath(genpath(spmPath));

    cd(InpaintingPath);
    subjects = dir('s*');

    for i = 1:length(subjects)
        fHeader = spm_vol(fullfile(HCPDatapath,subjects(i).name,'T1w','new_brainmask.nii'));
        mask = spm_read_vols(fHeader);
        indices_wb = find(mask);
        save(fullfile(InpaintingPath,subjects(i).name,'fHeader.mat'),'fHeader');
        %save(fullfile(InpaintingPath,subjects(i).name,'indices_wb.mat'),'indices_wb');
    end
    cd('/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine/connectivity-pipeline');
end