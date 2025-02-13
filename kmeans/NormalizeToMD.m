function NormalizeToMD(Path,prefix,deformfield)
% %     Inputs: 
% %     Path -- Path to data (folder)
% %     prefix -- prefix of files you want to normalize
% %     deformation -- deformation field (y*) (folder + filename)
% % 
% %     Normalize the volume to MNI space
    
    disp('Normalizing results to MNI...');
    fVolsFNlist = dir(fullfile(Path,[prefix,'*']));
    fVolsFNlist = struct2cell(fVolsFNlist); 
    fVolsFNlist = fVolsFNlist(1,:);
    fVolsFNlist_fullpath=cellfun(@(x) fullfile(Path,x), fVolsFNlist,'UniformOutput',false);
    
    spm_jobman('initcfg'); 
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deformfield};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = fVolsFNlist_fullpath';
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN
                                                            NaN NaN NaN];                                                  
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'MNI_';
    spm_jobman('run',matlabbatch);
end
    