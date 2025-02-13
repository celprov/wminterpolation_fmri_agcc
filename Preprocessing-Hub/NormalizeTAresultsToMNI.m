%% Normalize subjects to MNI before feeding them to TA+iCAPs pipeline 
%       for the extraction of the iCAPs

clear all; close all;

loc = '/home/atarun/fMRI-Sleep-Whole-Scan';
addpath(fullfile(loc, 'Preprocess', 'spm12'))

subjects = dir(fullfile(loc, 's*'));
TA= 'ta';
switch TA
    case 'ta'
        TAfolder = 'TotalActivation'
        
    case 'thresh'
        TAfolder = 'Thresholding/Alpha_1_99_Fraction_0DOT1'
end
files = {'Activity_inducing', 'Activity_related', 'Innovation'};
for i =[1,2,4:length(subjects)]  
    k=3
    overallfile = [];
    subject = subjects(i).name;
    display(strcat('running subject ...', subject))
    databasePath = fullfile(loc,subject,'TA_results');
    dirs = dir(fullfile(databasePath, strcat('Sleep_Data_',subject,'_subfolder_regressedOut_*')));
    structpath = fullfile(loc, subject, 'struct','Segmented');
    for j = 1:length(dirs)
        folder = fullfile(databasePath, dirs(j).name,TAfolder);
        newfolder = fullfile(folder, strcat(files(k), '_3D'));
        if k~=3
            TAfiles = dir(fullfile(newfolder{1,1},cell2mat(strcat(files(k),'_*'))));
        else
            TAfiles = dir(fullfile(newfolder{1,1}, 'Innovation_*'));
        end 
        l = length(TAfiles);  
       
           if l>100
                overlap = 10;
                if j==1
                    range = 1:l-overlap;
                elseif j==length(dirs)
                    range = overlap+2:l;
                else
                    range = overlap+2:l-overlap;
                end
                for m = 1:length(range)
                    filepath = fullfile(newfolder{1,1},TAfiles(range(m)).name);
                    overallfile = cat(1,overallfile,filepath);
                end
           end
        
    end
    overall = struct([]);
    for p=1:size(overallfile,1)
        overall{p} = overallfile(p,:);
    end

%     

    
    %% Normalize data to MNI
    spm_jobman('initcfg'); 
    deformfield = dir(fullfile(structpath, 'y_*'))
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {fullfile(loc,subject,'struct','Segmented',deformfield(1).name)};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = overall'
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];                                                  
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('run',matlabbatch);

    
    
    %%   Concatenate Normalized 
    clear overallfile
    clear filepath
    overallfile = [];
    Allvolumes = [];
    for j = 1:length(dirs)
        folder = fullfile(databasePath, dirs(j).name,TAfolder);
        newfolder = fullfile(folder, strcat(files(k), '_3D'));
        if k~=3
            TAfiles = dir(fullfile(newfolder{1,1},cell2mat(strcat('w',files(k),'_*'))));
        else
            TAfiles = dir(fullfile(newfolder{1,1}, 'wInnovation_*'));
        end 
%         if l>100
%                 overlap = 10;
%                 if j==1
%                     range = 1:l-overlap;
%                 elseif j==length(dirs)
%                     range = overlap+2:l;
%                 else
%                     range = overlap+2:l-overlap;
%                 end
%                 for m = 1:length(range)
%                     filepath = fullfile(newfolder{1,1},TAfiles(range(m)).name);
%                     overallfile = cat(1,overallfile,filepath);
%                 end

        
        for m = 1:length(TAfiles)
        filepath = fullfile(newfolder{1,1},TAfiles(m).name);
        
        overallfile = cat(1,overallfile,filepath);
        end
        
   
     
    end
    overall = struct([]);
    for p=1:size(overallfile,1)
        overall{p} = overallfile(p,:);
        Avol = spm_vol(overall{p});
        Avol2 = spm_read_vols(Avol);
        Allvolumes(:,:,:,p)=Avol2;
    end
    
    
     wInnovation = Allvolumes;
     whereto = fullfile(databasePath, 'Concatenated_Sleep_Data_TAresults_regressedOut',...
        'Thresholding')
    
    save(fullfile(whereto, 'wInnovation.mat'),'wInnovation', '-v7.3')
%      Convert to 4D and save
%     name_file = fullfile(databasePath, strcat('w',files(k), '.nii'));
%     spm_jobman('initcfg');
%     matlabbatch{1}.spm.util.cat.vols = overall';
%     matlabbatch{1}.spm.util.cat.name = 'wInnovation.nii';%char(strcat(files(k),'.nii'));%name_file{1,1};
%     matlabbatch{1}.spm.util.cat.dtype = 0;
%     spm_jobman('run',matlabbatch);
%     
%     from = fullfile(databasePath, dirs(1).name, 'TotalActivation',strcat(files(k),'_3D'),...
%         strcat('w',files(k),'.nii'));
%     movefile(cell2mat(from), cell2mat(fullfile(whereto,strcat('w', files(k),'.nii'))))
end
