
 % Function to normalize functional volumes to MNI
 % A.T. (02.18.2017)


%% Inputs:
%     funcPath - path of the functional volumes to normalize
%     segPath - path of the deformation field
%     prefix  - prefix of the volumes to be processed
%     MNIsize - size of the output voxel
 
function [] = NormalizeToMNI_size(funcPath,prefix,deformfield,MNIsize)

    fVolsFNlist = dir(fullfile(funcPath,[prefix,'*']));
    fVolsFNlist = struct2cell(fVolsFNlist); 
    fVolsFNlist = fVolsFNlist(1,:);
    fVolsFNlist_fullpath=cellfun(@(x) fullfile(funcPath,x), fVolsFNlist,'UniformOutput',false);

    
    spm_jobman('initcfg');
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deformfield};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = fVolsFNlist_fullpath';
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN
                                                            NaN NaN NaN];                                                  
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [MNIsize MNIsize MNIsize];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = strcat('w',int2str(MNIsize));
    spm_jobman('run',matlabbatch);
    
end
