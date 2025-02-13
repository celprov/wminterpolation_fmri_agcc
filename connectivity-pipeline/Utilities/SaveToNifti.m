function [] = SaveToNifti(param, Utr, cases)

   
        wheretosave = fullfile(param.SaveDirectory, cases);
        if ~exist(wheretosave, 'dir')
            mkdir(wheretosave)
        end

        fHeader = spm_vol(fullfile(param.HCPDatapath,param.subject, 'T1w','extracted','new_brainmask.nii'));
        hdr=cbiReadNiftiHeader(fHeader.fname);

        
        indices = param.G.indices_wb;
       
        for i = 1:size(Utr,2)

            V = zeros(fHeader.dim);
            V(indices) = Utr(:,i);
            niftiFile = fullfile(wheretosave,[sprintf('Eigenmode_%.4d',i),'.nii']);

            cbiWriteNifti(niftiFile,V,hdr,'float32');
        end
    
        
end